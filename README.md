pubmedBatch in Python (IRDC branch)
==============================
Send your searches in batch against Pubmed!

Installation
------------
You will need `python`, `flask`, `mongodb` to make the app work

And to speed up the app, you should download the following dtd files
```bash
wget -O ~/.biopython/Bio/Entrez/DTDs/nlmmedlinecitationset_160101.dtd http://www.ncbi.nlm.nih.gov/corehtml/query/DTD/nlmmedlinecitationset_160101.dtd
wget -O ~/.biopython/Bio/Entrez/DTDs/pubmed_160101.dtd http://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_160101.dtd
wget -O ~/.biopython/Bio/Entrez/DTDs/bookdoc_160101.dtd http://www.ncbi.nlm.nih.gov/corehtml/query/DTD/bookdoc_160101.dtd
wget -O ~/.biopython/Bio/Entrez/DTDs/efetch.dtd http://eutils.ncbi.nlm.nih.gov/eutils/dtd/20131226/efetch.dtd
```
Configuration
-------------
The configuration is set in `config.cfg`.
You can change the settings according to your use.
By default, the app is running on `localhost:8000/pubmedbatch`

Introduction
------------
if you have a list of genes to search against Pubmed, you can put them into a column together with other relevant information, and upload it as a `csv` file. **Note** that the `csv` file requires a header.

The app will try to calculate a `Pubmed score` based on the occurrence of the OR terms in the titles/abstracts, and sort articles / genes based on the *score per article* and *score per gene*, respectively.

It recognises some column headers, such as 'Func', 'ExonicFunc', 'Pred' in order to heuristically aggregate the predicting damaging effect of a given variant:  
* If `splic` | `stop` | `frame` | `del` | `insert` found, set Pred score to 1000
* Else, each `D` or `A` gives 10, `P` gives `5`, `C` gives 6, `T` | `B` | `N` gives -1  

It also recognises column with header 'FILTER'. If 'FILTER' is `PASS`, the `del row` button is coloured green. If it is `FAIL`, the `del row` button is coloured yellow. Else the button is coloured red.

Features
--------
* It saves the results of each search in a database. It won't search the same term in Pubmed again until 30 days later (can be changed in `config.cfg`). After 30 days, it will only update the search results by searching Pubmed for articles published in the recent 30 days.
* The table is sortable by any column.
* The Pubmed ID is clickable.
* You can provide a gene list (space separated) to highlight genes.
* You can provide a blacklist to skip genes.
* You can use the files in the *testfile* folder to play with the app.
* One caveat of the app is, deleting row won't delete it in the actual data. 

IRDC features
-------------
* It refers to [RetNet](https://sph.uth.edu/Retnet/) to mark genes as **D** for *dominant*, **R** for *recessive*, **X** for *X-linked*, **M** for *mitochondrion* and **U** for *unknown*.
* Click on the marked genes will pop disease details and available OMIM links.
* In the `scraper` folder, `scrape_RetNet.py` can scrape the content of [RetNet](https://sph.uth.edu/Retnet/disease.htm) and store it in JSON format.

Future work
-----------
* Will take tab delimited files