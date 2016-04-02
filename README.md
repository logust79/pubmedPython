pubmedBatch in Python
==============================
Send your searches in batch against Pubmed!

Installation
------------
You will need python and flask to make the app work

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
