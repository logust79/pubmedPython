#!/usr/bin/env python
from Bio import Entrez
import itertools
import json
import os
import pymongo
import logging
import sys
from utils import *
from urllib2 import HTTPError
from flask import Flask, request, session, g, redirect, url_for, abort
from flask import render_template, flash, jsonify, Response, send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
from flask_errormail import mail_on_500
from flask import Response
from collections import defaultdict, Counter
from werkzeug.contrib.cache import SimpleCache
from multiprocessing import Process
import glob
import sqlite3
import traceback
import time
from flask_debugtoolbar import DebugToolbarExtension
from werkzeug.exceptions import default_exceptions, HTTPException
import csv
from functools import wraps
import StringIO
from urlparse import urlparse
import pickle
import re
import xml.etree.ElementTree as xml

from flask.ext.session import Session

logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

ADMINISTRATORS = ('jing.yu@ndcn.ox.ac.uk',)
app = Flask(__name__)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True

# Load default config and override config from an environment variable
app.config.from_pyfile('config.cfg')

# Check Configuration section for more details
SESSION_TYPE = 'mongodb'
app.config.from_object(__name__)
sess=Session()
sess.init_app(app)

# read RetNet file
retnet_f = 'retnet.json'
RETNET = json.load(open(retnet_f, 'r'))

# progressbar
'''
{
    'random_p_id':{
        'total':456,
        'count':123,
        'status':['running','done']
    },
    ...
}
'''
PROGRESS_BAR = {}

'''
initiate a progress instance
arg: total length of genes
return: progress_id
'''

def init_progress_bar(id,length):
    # check id
    if id in PROGRESS_BAR:
        if PROGRESS_BAR[id]['status'] != 'done':
            return 'the id already exists in PROGRESS_BAR'

    # initialise progress_bar
    PROGRESS_BAR[id] = {
            'total': length,
            'count':0,
            'message': '',
            'status':'running'
    }
    return 0

'''
update progress
arg: {
    id: id, 
    message: message,
    step: 1
    }
default step 1
'''

def update_progress_bar(obj):
    # check if id in PROGRESS_BAR
    if not obj['id'] in PROGRESS_BAR:
        return 'ID does not exist in PROGRESS_BAR'

    # update progress
    if not 'step' in obj:
        obj['step'] = 1
    PROGRESS_BAR[obj['id']]['count'] += obj['step']

    PROGRESS_BAR[obj['id']]['message'] = obj['message']
    # done?
    if PROGRESS_BAR[obj['id']]['count'] == PROGRESS_BAR[obj['id']]['total']:
        PROGRESS_BAR[obj['id']]['status'] = 'done'

'''
kill a progress
'''

def kill_progress_bar(key):
    if key in PROGRESS_BAR:
        del PROGRESS_BAR[key]

'''
to check if an iterable is empty
'''
def peek(iterable):
    try:
        first = next(iterable)
    except RuntimeError:
        return None
    except StopIteration:
        return None
    return first, itertools.chain([first], iterable)


'''
find the freaking PID, Title or Abstract no matter what!
'''
def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item

def connect_db(dbname=None):
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    if not dbname: dbname=app.config['DB_NAME']
    print(dbname)
    return client[dbname]

def get_db(dbname=None):
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if dbname is None: dbname=app.config['DB_NAME']
    if not hasattr(g, 'db_conn'):
        g.db_conn=dict()
        g.db_conn[dbname] = connect_db(dbname)
    elif dbname not in g.db_conn:
        g.db_conn[dbname] = connect_db(dbname)
    return g.db_conn[dbname]

"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def scrutinise(obj):
    if obj['lag']:
        obj['lag'] = obj['lag']/3600/24 # convert it to days
            # need to update
        search_results = Entrez.read(Entrez.esearch(db='pubmed', term=obj['smashed_all'], reldate=obj['lag'], datetype='pdat', usehistory='y'))
    else:
        # just search
        search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=obj['smashed_all'], usehistory='y'))
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'total_score':0}
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed",
                                   restart=0,
                                   retmax=50,
                                   retmode="xml",
                                   webenv=search_results['WebEnv'],
                                   query_key=search_results['QueryKey']
                                   )
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print('Received error from server %s' % err)
            else:
                print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    if peek(record):
        # got something. let's do some calculation
        for r in record:
            # calculate score
            score = 0
            pid = str(find_item(r, 'PMID'))
            abstract_list = find_item(r, 'AbstractText')
            # parse abstract
            abstract = ''
            if abstract_list:
                for a in abstract_list:
                    if hasattr(a, 'attributes') and 'Label' in a.attributes:
                        abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                        abstract = abstract + a + '<br/>'
                    else:
                        abstract = abstract + a

            title = find_item(r, 'ArticleTitle')
            if title:
                score = score + len(obj['reg'].findall(title))
            if abstract:
                score = score + len(obj['reg'].findall(abstract))

            # add result to genes[gene_name]
            if score:
                results['results'].append({
                    'id': pid,
                    'title': title,
                    'abstract': abstract,
                    'score': score
                })
                results['total_score'] = results['total_score'] + score
    results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results

def get_pred_score(obj):
    # for the batch_pubmed route.
    # calculate the pred score
    # [D/A].each = 10, [P].each = 5, [C].each = 6, [T/B/N].each = -1. If there is a splicing/insertion/deletion event, the score is set as 1000. Not given is set as 0
    # ref: https://github.com/plagnollab/DNASeq_pipeline/blob/master/GATK_v2/filtering.md
    pred = 0
    if ('Func' in obj and re.search('splic', obj['Func'])) or ('ExonicFunc' in obj and re.search(r'stop|frame|del|insert', obj['ExonicFunc'])):
        pred = 1000;
    else:
        for k in obj:
            if re.search('Pred', k):
                if obj[k] == 'D' or obj[k] == 'A':
                    pred = pred + 10
                elif obj[k] == 'P':
                    pred = pred + 5
                elif obj[k] == 'C':
                    pred = pred + 6
                elif obj[k] == 'T' or obj[k] == 'B' or obj[k] == 'N':
                    pred = pred - 1
                else:
                    pass
    return pred;

@app.route('/pubmedbatch/', methods=['GET', 'POST'])
#@requires_auth
def pubmedbatch_main():
    # this is the main page of pubmedBatch
    # It allows a user to create a "folder" to hold results in the mangodb. Post request will handle the 'folder' creation
    '''
    db.results = {
            'user_id':'Jing',
            'folder':{
                foo':{
                    'test':{JSON}
                    'test(1):{JSON}
            },
        },
    }
    '''
    user = session.get('user') or app.config['DEFAULT_USER']

    db = get_db('pubmedbatch')
    user_db = db.results.find_one({'user_id':user},{'_id':0})

    if request.method == 'POST':
        # time to create a folder
        # first check if the folder already exists. Yes? pass. No? create
        folder = request.form['create-folder']

        user_folders =  [d for d in user_db['folder']]
        if folder not in user_folders:
            # set dot.notation
            user_folder = 'folder.' + folder
            db.results.update({'user_id':user}, {'$set':{user_folder:{} }})
    else:
        # get folders. If user not exists in pubmed db, create it
        #print(res)
        if not user_db:
            # A new user dropping by...
            print('A "user" is being made in pubmedDB!')
            db.results.insert({'user_id':user, 'folder':{}})

    # let's render the template  
    user_db = db.results.find_one({'user_id':user},{'_id':0})

    return render_template( 'pubmedbatch_main.html', 
            user = user,
            folders = [d for d in user_db['folder']]
    )

'''
delete a folder in pubmedbatch
'''
@app.route('/pubmedbatch_folderdel/<folder>', methods=['POST'])
#@requires_auth
def pubmedbatch_delfolder(folder):
    user = session.get('user') or app.config['DEFAULT_USER']

    db = get_db('pubmedbatch')
    user_db = db.results.find_one({'user_id':user},{'_id':0})
    folders = [d for d in user_db['folder']]
    if folder in folders:
        folder = 'folder.' + folder
        db.results.update({'user_id':user}, {'$unset': {folder:''}})
        return redirect('/pubmedbatch/')
    else:
        return 'Folder: <b>%s</b> does not exist!' % folder

'''
main business
'''
@app.route('/pubmedbatch/<folder>', methods=['GET', 'POST'])
#@requires_auth
def pubmedbatch(folder):
    # This is the main business
    user = session.get('user') or app.config['DEFAULT_USER']
    db = get_db('pubmedbatch')
    if request.method == 'POST':
        # post. get form data, return JSON
        ##########################################################
        # get form data
        #########################################################
        column = int(request.form['column']) - 1
        Entrez.email = request.form['email']
        AND_term = request.form['AND']
        OR_term = request.form['OR']
        verbose = request.form.get('verbose','')

        # set session
        session['AND'] = AND_term
        session['OR'] = OR_term
        session['EMAIL'] = request.form['email']
        session['COLUMN'] = column + 1
        
        csv_file = request.files['csv_upload']
        file_name = csv_file.filename
        known_genes = request.files['known_genes'] or ''
        mask_genes = request.files['mask_genes'] or ''
        #########################################################
        # parse known genes and mask genes
        known_genes = known_genes.read().split() if known_genes else ''
        mask_genes = mask_genes.read().split() if mask_genes else ''
        #########################################################
        # read csv. has to make the csv string looking like a filehandle
        csv_content = csv_file.read()

        '''
        progress bar initiation
        '''
        # how many lines? minus header
        csv_lines = len(re.findall(r'\n', csv_content)) - 1
        # a random id
        progress_id = user + request.form['rid']
        bad = init_progress_bar(progress_id, csv_lines)
        if bad:
            print 'weird, progress id conflict!'
        '''
        '''
        csvreader = csv.reader(StringIO.StringIO(csv_content), delimiter=',', quotechar='"')
        # life time from config
        life = app.config['PUBMEDBATCH_LIFE']
        # number of lines?
        line_num = len(re.findall('\n', csv_content))
        # format terms
        OR = OR_term.split()
        OR.sort()
        # make reg term for later use (calculating pubmed score
        reg = '\\b|\\b'.join(OR)
        reg = '\\b' + reg + '\\b'
        reg = re.compile(reg, re.IGNORECASE)

        AND = AND_term.split()
        AND.sort()
        smashed_OR = ['"' + t + '"' + '[Title/Abstract]' for t in OR]
        smashed_AND = ['"' + t + '"' + '[Title/Abstract]' for t in AND]
        smashed_OR = ' OR '.join(smashed_OR)
        smashed_AND = ' AND '.join(smashed_AND)
        smashed_term = ' AND (' + smashed_OR + ')'
        if smashed_AND:
            smashed_term += ' AND ' + smashed_AND
        ###########################################################
        # it's time to read the csv file
        ###########################################################
        row_num = -1
        header = [] # header
        output = [] # all the results get pushed here
        for row in csvreader:
            row_num += 1
            # read header
            if not row_num:
                header = row
                # replace . with _ in header
                header = [re.sub('\.','&#46;', h) for h in header]
                # add 2 columns after HUGO
                header[column+1:column+1] = ['ref(pubmedID)', 'pubmed_score']
                # add a pred score at the beginning
                header[0:0] = ['pred_score']
                continue
            # read in real data
            gene_name = row[column]
            if gene_name == 'NA':
                continue
            # get rid of any parentheses and their content
            gene_name = re.sub(r'\([^)]*\)?','',gene_name)

            # update progress bar
            update_progress_bar({'id':progress_id, 'message': '%s;' % gene_name})

            genes={} # storing 
            print gene_name
            smashed_all = gene_name + smashed_term
            ####################################################
            # first to see if masked, pass
            # else
            #   db.cache
            #       when in db.cache and up-to-date, get the pubmed result
            #       when in db.cache and out-of-date, add new pubmed result
            #       when not in db.cache, search pubmed and store
            ####################################################
            # db.cache's structure:
            #  {['key': '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()]), 
            #   'result': {pubmed_result},
            #   'date': now]}
            if gene_name in mask_genes:
                # masked. don't search, just return the row
                if not verbose:
                    continue
                row[column+1:column+1] = [{'total_score': 0, 'results': ['masked']}, 0]
                # add a placeholder for pred_score
                row[0:0] = [0]
                ha = {}
                for col in range(len(header)):
                    ha[header[col]] = row[col]
                ha['pred_score'] = get_pred_score(ha)
                output.append(ha)
            else:
                # not masked
                now = time.mktime(time.localtime()) #get time in seconds
                lag = 0 # use it as a flag of how to search. 0 = search; now-saved['date'] = update; 
                term = '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()])
                # now check if the result in the db is uptodate
                saved = db.cache.find_one({'key': term})
                if saved:
                    lag = now - saved['date']
                    # has record. let's see if it is out of date
                    if lag  <= life:
                        # not out of date. let's use it
                        genes[gene_name] = saved['result']

                if gene_name not in genes:
                    # need to search
                    # handle = Entrez.einfo()
                    # record = Entrez.read(handle)
                    # handle.close()
                    results = scrutinise({'lag':lag, 'smashed_all': smashed_all, 'reg':reg})

                    # update the database, and maybe results
                    if lag:
                        # need to update the database
                        # old result = saved
                        results['total_score'] = results['total_score'] + saved['result']['total_score']
                        results['results'].extend(saved['result']['results'])
                        # update database now
                        db.cache.update({'key': term}, {'$set': {'result': results, 'date': now}})
                    else:
                        # write to the database
                        db.cache.insert({'key': term, 'result': results, 'date': now})

                    genes[gene_name] = results

                # verbose and no score?
                if not (verbose or genes[gene_name]['total_score']):
                    continue
                else:
                    # known genes?
                    if gene_name in known_genes:
                        genes[gene_name]['known'] = 1
                        # give the rest a minimal score to keep them on the list
                        genes[gene_name]['total_score'] = max(1, genes[gene_name]['total_score'])
                    else:
                        genes[gene_name]['known'] = 0
                    # Retnet?
                    if gene_name in RETNET:
                        genes[gene_name]['disease'] = RETNET[gene_name]['disease']
                        genes[gene_name]['omim'] = RETNET[gene_name]['omim']
                        genes[gene_name]['mode'] = RETNET[gene_name]['mode']
                        # reassign total score according to mode
                        if RETNET[gene_name]['mode'] == 'd' or RETNET[gene_name]['mode'] == 'x':
                            genes[gene_name]['total_score'] = max(100, genes[gene_name]['total_score'])
                        elif re.search(r'/(?!dominant)/', AND_term):
                            # not searching doimant, also assign others to 100
                            genes[gene_name]['total_score'] = max(100, genes[gene_name]['total_score'])
                        else:
                            # give the rest a minimal score to keep them on the list
                            genes[gene_name]['total_score'] = max(1, genes[gene_name]['total_score'])
                    
                    # add pubmed result
                    row[column+1:column+1] = [genes[gene_name], genes[gene_name]['total_score']]
                    # add a placeholder for pred_score
                    row[0:0] = [0]
                    ha = {}
                    for col in range(len(header)):
                        ha[header[col]] = row[col]
                    ha['pred_score'] = get_pred_score(ha)
                    output.append(ha)

        # save results to folder
        # if the file name already exists?
        #folders = db.results.find_one({'user_id':user},{'_id':0})['folder']
        the_folder = db.results.find_one({
            'user_id': user},{'_id': 0})['folder'][folder]
        files = [ f for f in the_folder]
        file_name = os.path.splitext(file_name)[0]
        num = 0
        while file_name in files:
            num += 1
            file_name = re.sub(r'(\(\d+\))?$', '(%s)' % num, file_name)
    
        # get the search term, to display
        search_term = 'AND[ <b>%s</b> ]; OR[ <b>%s</b> ]' % (', '.join(AND), ', '.join(OR))
        # write to database and win
        final_result = [header, output, search_term]
        the_file = 'folder.' + folder + '.' + file_name
        db.results.update({'user_id': user}, {
            "$set": {the_file: final_result}
        })
        return json.dumps([header, output, search_term, file_name])
    else:
        # get. display page
        # First see if folder exists. if not, return error
        user_folders = db.results.find_one({
            'user_id': user},{'_id': 0})['folder']
        # get AND an OR field
        AND = session.get('AND') or ''
        OR = session.get('OR') or app.config['PUBMEDBATCH_OR']
        EMAIL = session.get('EMAIL') or app.config['PUBMED_EMAIL']
        COLUMN = session.get('COLUMN') or '4'
        #user_folders = user_folders['folder']
        if folder not in user_folders:
            return "Error: " + folder + " does not exist!"
        # get the files in the folder to display
        files = [e for e in user_folders[folder]]
        return render_template('pubmedbatch.html',
            home_pubmedbatch = home_pubmedbatch,
            files = files,
            user = user, 
            folder = folder,
            AND = AND,
            OR = OR,
            EMAIL = EMAIL,
            COLUMN = COLUMN
        )


@app.route('/pubmedbatch_rename/<folder>', methods=['POST'])
# @requires_auth
def pubmedbatch_rename(folder):
    '''
    Rename a file in pubmedbatch
    '''
    user = session.get('user') or app.config['DEFAULT_USER']

    new_name = request.form['new-name']
    old_name = request.form['old-name']

    # check if new-name already exists
    db = get_db('pubmedbatch')
    the_folder = db.results.find_one({
            'user_id': user}, {'_id': 0})['folder'][folder]
    files = [f for f in the_folder]
    # sanity check
    if new_name in files:
        # throw an error
        return '''<span style="color:orange">Oops, the new name:
             <b>%s</b> already exists in the folder</span>''' % new_name
    if old_name not in files:
        # shouldn't happen though
        return '''<span style="color:orange">Oops, the old name: <b>%s</b>
        does not exist in the folder, somehow</span>''' % old_name
    # rename
    # mongodb doesn't like multiple '$' for the time being,
    #  so need to find hard index
    old_name = 'folder.' + folder + '.' + old_name
    new_name = 'folder.' + folder + '.' + new_name
    db.results.update({'user_id': user}, {'$rename': {old_name: new_name}})

    return 'The file has been successfully renamed!'


@app.route('/pubmedbatch_del/<path:path>', methods=['POST'])
# @requires_auth
def pubmedbatch_del(path):
    '''
    Delete a file in pubmedbatch
    '''
    user = session.get('user') or app.config['DEFAULT_USER']

    # get folder and file from route
    m = re.search(r'([^/]+)/([^/]+)', path)
    (folder_name, file_name) = m.groups()
    file_name = re.sub('%20', ' ', file_name)

    db = get_db('pubmedbatch')
    db.results.update({'user_id': user}, {'$unset': {
                      'folder.%s.%s' % (folder_name, file_name): ''}})

    return '1'


@app.route('/pubmedbatch/<path:path>', methods=['POST'])
# @requires_auth
def pubmedbatch_getfile(path):
    '''
    get a file in pubmedbatch
    '''
    user = session.get('user') or app.config['DEFAULT_USER']

    # get folder and file from route
    m = re.search(r'([^/]+)/([^/]+)', path)
    (folder_name, file_name) = m.groups()
    file_name = re.sub('%20', ' ', file_name)
    # get file from database
    db = get_db('pubmedbatch')
    file = db.results.find_one({'user_id':
            user})['folder'][folder_name][file_name]

    return json.dumps(file)


@app.route('/pubmedbatch_progress_bar/<id>')
def pubmedbatch_progress(id):
    '''
    progress bar query
    '''
    user = session.get('user') or app.config['DEFAULT_USER']
    progress_id = user + id
    return json.dumps(PROGRESS_BAR[progress_id])


if __name__ == "__main__":
    home = ''
    home_pubmedbatch = '/pubmedbatch'
    app.run(host='0.0.0.0', port=8000, threaded=True)
