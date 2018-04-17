# pubmed-articles-grabber
Grabs and downloads full-text versions of PubMed records in different formats.

Grabs and downloads full-text versions of PubMed records in different formats. It arranges them in folders named after the file formats (XML, PDF, ..). It is recommended to have two extra files: one containing the click through token to access articles that needs subscription, and the other containing the list of wanted PMIDs.  The names of the files are 'clickThroughToken.txt' and 'wanted.csv', respectively.  

The contents of:      
- wanted.csv: rows of PMIDs with 'PMID' as a header(i.e. first row).     
- clickThroughToken.txt: the researcher's Click-through-token, obtained as in https://support.crossref.org/hc/en-us/articles/115001376963-Click-through-service-for-researchers. 

The main functions to grab the articles  are: `grab()`, `grabViaCrossref()` and `grabViaPMCOAI()`. The rest of the functions are seervice functions that can be used for putposes other than the main purpose(e.g. converting XML to txt formats).    
    
###     Recommended steps::  
1. Put the list of PMIDs you want to "grab" in the 'wanted.csv' file in a row by row basis.    
2. (optional) Put your click-through-token in the 'clickThroughToken.txt' file. This step is optional but some publishers require having a click-through-token.    wanted
3. Create an object of this class and use its main grabber functions. Below is an example of its use (also can be found in 'main.py').    
4. Do not forget to change the email address from "example@example.domain" to yours.    
    
##  Usage::        
```
from pubMedArticleGrabber import PubMedArticleGrabber    
wanted = PubMedArticleGrabber('wanted', 'example@example.domain')
wanted.grab()  
```
