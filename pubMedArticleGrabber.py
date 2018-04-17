try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import os
import sys
import time
import pickle
import re
import csv
import unicodecsv
import pandas as pd
import requests
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from ssl import CertificateError
from collections import defaultdict
from Bio import Entrez
from habanero import Crossref

class PubMedArticleGrabber(object):
    """
    Grabs and downloads full-text versions of PubMed records in different formats. It arranges them in folders named after the file formats (XML, PDF, ..). It is recommended to have two extra files: one containing the click through token to access articles that needs subscription, and the other containing the list of wanted PMIDs.  The names of the files are 'clickThroughToken.txt' and 'wanted.csv', respectively.

    The contents of: 
    - wanted.csv: rows of PMIDs with 'PMID' as a header(i.e. first row).
    - clickThroughToken.txt: the researcher's Click-through-token, obtained as in https://support.crossref.org/hc/en-us/articles/115001376963-Click-through-service-for-researchers 

    The main functions to grab the articles  are: grab, grabViaCrossref and grabViaPMCOAI. The rest of the functions are seervice functions that can be used for putposes other than the main purpose(e.g. converting XML to txt formats).

    Usage::
    from pubMedArticleGrabber import PubMedArticleGrabber
    
    wanted = PubMedArticleGrabber('wanted', 'xxx@xxx.xx')
    wanted.grab()
"""

    def __init__(self, pmidsListOrCSVfile, email):
        """
        Constructor. Used for initialization.

        :param pmidsListOrCSVfile: [list, set or string] a list (or a set) of PMIDs, or a name of a CSV file  containg  the list of PMIDs with the header 'PMID'.
        :param email: [string] an email address for Entrez and Crossref in case of need to contact you (e.g. send warnigs about download limits).

        Usage::
        from pubMedArticleGrabber import PubMedArticleGrabber
        wanted = PubMedArticleGrabber('name of csv file', 'xxx@xxx.xx')
        others = PubMedArticleGrabber(['1047458', '1050021'], 'xxx@xxx.xx')
        someOthers = PubMedArticleGrabber({'1047458', '1047458', '1050021'}, 'xxx@xxx.xx')
        """
        self.pmids = pmidsListOrCSVfile
        Entrez.email = email
        cr = Crossref(mailto=email)


    def grab(self):
        """The main function to grab articles via both PubMedCentral and Crossref"""
        self.grabViaPMCOAI()
        self.grabViaCrossref()

    def grabViaCrossref(self):
        """Grabs articles via Crossref, after converting PMIDs to DOIs through Entrez"""
        self.pmids

        # converting pmids to dois to use them for downloading full text articles
        self.pmidTodoi(self.pmids,'dois')

        # reading the DOIs
        dois = pd.read_csv('dois.csv')
        dois = set(dois['DOI'])

        # the variable "articles" is a list of dictionaries (a dictionary for each article). The dictionary is the result of querying habanero crossref using the article's doi. Here we are trying to load the variable "articles" if it was already pickled or start with an empty list.
        try:
            with open('articles.pkl', 'rb') as f:
                articles = pickle.load(f)
        except FileNotFoundError:
            articles = list()

        # some dois are not found in crossref. Here we are trying to read them or start with an empty set
        try:
            notFounddois = set(self.csvTolist('full text/doi/dois notFoundin crossRef'))
        except FileNotFoundError:
            notFounddois = set()

        # querying crossref by the dois and saving the results
        # we wrap the main process in a while loop to resume the process if it was interrupted due to network errors.
        firstTime = True
        while(dois):
            # first we filter out the dois that were already processed
            collecteddois = set()
            for article in articles:
                collecteddois.add(article['message']['DOI'])
            collecteddois = {self.extractDOI(i).lower() for i in collecteddois}
            notFounddois = {self.extractDOI(i).lower() for i in notFounddois}

            dois = {i.lower() for i in dois}
            dois = dois - notFounddois - collecteddois
            if firstTime:
                lenOfIn = len(dois)
                firstTime = False
            elif len(dois) == lenOfIn:
                print('Remaining: ', len(dois) ,'DOIs')
                break

            # This is the main process of querying crossref for the dois that are not processed yet
            for doi in dois:
                try:
                    articles.append(cr.works(ids=doi))
                except (HTTPError, requests.exceptions.HTTPError, URLError, CertificateError, requests.exceptions.SSLError, TimeoutError) as e:
                    if 'Not Found' in str(e):
                        notFounddois.add(doi)
                    else:
                        print(doi)
                        raise e
        notFounddois = {self.extractDOI(i).lower() for i in notFounddois}
        self.listTocsvRows(notFounddois, 'full text/doi/dois notFoundin crossRef')

        # remove duplicate results
        temparticles=list()
        seendois = set()
        for article in articles:
            if article['message']['DOI'] not in seendois:
                temparticles.append(article)
                seendois.add(article['message']['DOI'])
        articles = temparticles
        del temparticles

        # pickle the variable "articles", due to its size and to the time (network time) needed to build it.
        with open('articles.pkl', 'wb') as f:
            pickle.dump(articles, f)
        with open('articles.pkl', 'rb') as f:
            articles = pickle.load(f)

        # Extracting urls that enable us to download full text articles. Some dois does not have a url in crossref
        articlesURLs = dict()
        doisWithNoLink = set()
        for article in articles:
            doi = article['message']['DOI']
            try:
                links = article['message']['link']
            except KeyError:
                doisWithNoLink.add(doi)
                continue
            urls = set()
            [urls.add((link['URL'], link['content-type'])) for link in links]
            if urls:
                articlesURLs[doi] = urls
            else:
                doisWithNoLink.add(doi)
        self.listTocsvRows(doisWithNoLink, 'full text/doi/dois With No Link')


        # Read dois and urls that are already downloaded if they exist. Also read urls that are known to be broken
        try:
            downloadabledois = set(pd.read_csv('full text/doi/downloadable dois.csv')['DOI'].dropna())
        except FileNotFoundError:
            downloadabledois = set()
        try:
            downloadableURLs = set(pd.read_csv('full text/doi/downloadable dois.csv')['URL'].dropna())
        except FileNotFoundError:
            downloadableURLs = set()
        try:
            badURLs = set(pd.read_csv('full text/doi/bad urls.csv')['URL'].dropna())
        except FileNotFoundError:
            badURLs = set()

        # click-through-token is obtained from ORCID and is saved as a dictionary pickle. It is passed as a header with the url requests to acquire access to full text articles that requires subscription or existence of this token for data mining purposes.
        try:
            with open('clickThroughToken.txt', 'r') as f:
#                h = pickle.load(f)
                h = {'CR-Clickthrough-Client-Token': f.read()}
        except FileNotFoundError:
            print('There is no XRef click through token')
            h = None

        # downloading the articles
        for doi, urls_types in articlesURLs.items():
            # if the article was available from cambridge core
            if '10.1017/s' in doi:
                url = self.getCambridgeURL(doi)
                if url and url not in downloadableURLs|badURLs:
                    self.download(doi, url, h)
            else:
                for url_type in urls_types:
                    url, crContentType = url_type
                    if url and url not in downloadableURLs|badURLs:
                        # springer links of the form "springerlink.com/..." need to be fixed 
                        if 'springerlink' in url:
                            url = self.fixSpringerLinkURL(url)
                            if url and url in downloadableURLs|badURLs:
                                continue
                        self.download(doi, url, h)
        print('finish')

        #Retry (several times) downloading bad (non-downloaded) articles to avoid timeout and similar errors
        lenOfIn = 1
        lenOfOut = lenOfIn - 1
        while lenOfOut<lenOfIn:
            try:
                badArticles = pd.read_csv('full text/doi/bad urls.csv')
                badArticles = badArticles[['DOI', 'URL']].dropna()
                badURLs = set(badArticles['URL'].dropna())
            except FileNotFoundError:
                badArticles = pd.DataFrame(data={'DOI':[], 'URL':[]})
                badURLs = set()
            downloadableURLs = set(pd.read_csv('full text/doi/downloadable dois.csv')['URL'].dropna())
            lenOfIn = len(badURLs - downloadableURLs)
            # downloading the articles
            for index, row in badArticles.iterrows():
                doi = row['DOI']
                url = row['URL']
                if url == '###':
                    continue
                # if the article was available from cambridge core
                if '10.1017/s' in doi:
                    url = self.getCambridgeURL(doi)
                    if url and url not in downloadableURLs:
                        self.download(doi, url, h)
                elif url and url not in downloadableURLs:
                    # springer links of the form "springerlink.com/..." need to be fixed 
                    if 'springerlink' in url:
                        url = self.fixSpringerLinkURL(url)
                        if url and url in downloadableURLs:
                            continue
                    self.download(doi, url, h)
            print('finish')
            downloadableURLs = set(pd.read_csv('full text/doi/downloadable dois.csv')['URL'].dropna())
            lenOfOut = len(badURLs - downloadableURLs)
        print('Remaining: ', lenOfOut ,'IDs')


    def grabViaPMCOAI(self):
        """Grabs articles via PMC-OAI service, after converting PMIDs to PMCIDs through Entrez. It constructs a specific URL for each article to use the PMC-OAI service"""
        self.pmidTopmcid(self.pmids,'pmcids')
        d = self.csvTodict('pmcids')
        pmcids = {v.pop() for k, v in d.items()}
        pmcids.discard('PMCID')
        pmcids.discard('NoPMCID')
        noTxt = self.getFullText(pmcids)
        self.listTocsvRows(noTxt,'full text/txt/noTxt')

    def pmidTopmcid(self, pmidsListOrCSVfile, pmcidsCSVfile):
        # reading the pmids
        pmids = pmidsListOrCSVfile
        if type(pmids) is str:
            pmids = self.csvTolist(pmidsListOrCSVfile)
        try:
            if type(pmids[0]) is list and len(pmids[0]) == 2:
                pmids = [row[0] for row in pmids]
        except TypeError:
            pass
        # retreive pmcid for each pmid
        temp = pmids
        while temp.copy():
            for pmid in pmids.copy():
                handle = Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc", id=pmid)
                result = Entrez.read(handle)

                pmcid=''
                if len(result[0]['LinkSetDb']) == 1:
                    pmcid = result[0]['LinkSetDb'][0]['Link'][0]['Id']
                elif not result[0]['LinkSetDb']:
                    pmcid = 'NoPMCID'
                else:
                    pmcid = result[0]['LinkSetDb'][0]['Link'][0]['Id'][0]
                    pmcid += 'MoreThanPMCID'
                    print('MoreThanPMCID')
                    
                row = pmid, pmcid
                self.writecsvRow(pmcidsCSVfile, row, ['PMID', 'PMCID'])
                pmids.remove(pmid)
            temp = pmids


    def pmcidToXml (self, pmcid):
        url = 'https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:' +str(pmcid) + '&metadataPrefix=pmc'
        page = urlopen(url)
        xml = page.read()
        fileName = 'full text/pmcid/xml/' +str(pmcid)+ '.xml'
        directory = os.path.dirname(fileName)
        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(fileName, 'wb') as f:
            f.write(xml)


    def xmlToTxt (self, pmcid):
        noTxt = False
        fileName = 'full text/pmcid/xml/'+ str(pmcid)
        directory = os.path.dirname(fileName)
        if not os.path.exists(directory):
            os.makedirs(directory)

        tree = ET.parse(fileName + '.xml')
        root = tree.getroot()
        rootNS, _ = root.tag.split('}')
        rootNS = rootNS [1 :]
        nsmap = {'rootNS':rootNS}
        metadata = root.find('.//rootNS:metadata', nsmap)
        try:
            article = metadata[0]
        except TypeError:
            article = None
            noTxt = True
        if article:
            articleNS, temp = article.tag.split('}')
            assert temp == 'article', 'No article tag'
            del(temp)
            articleNS = articleNS[1 :]
            nsmap['articleNS'] = articleNS
            body = article.find('.//articleNS:body', nsmap)
            if body:
                txt = ' '.join(ET.tostringlist(body, encoding='unicode', method='text'))
                fileName = 'full text/pmcid/txt/' +str(pmcid)+ '.txt'
                directory = os.path.dirname(fileName)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                with open('full text/pmcid/txt/' +str(pmcid)+ '.txt', 'wb') as txtFile:
                    txtFile.write(txt.encode('utf-8'))
            else:
                noTxt = True

        return noTxt


    def getFullText(self, pmcids):
        noTxt=list()
        temp = pmcids
        while temp.copy():
            for pmcid in pmcids.copy():
                purePmcid = int(''.join(filter(str.isdigit, pmcid)))
                self.pmcidToXml(purePmcid)
                if self.xmlToTxt(purePmcid):
                    noTxt.append(purePmcid)
                    print(purePmcid)
                pmcids.remove(pmcid)
            temp = pmcids
        return noTxt


    def pmidEntrezSummaryRecordTodoi(self, esummaryRecord):
        summaryRecord = esummaryRecord[0]
        pmid = summaryRecord['Id']
        try:
            doi = summaryRecord['DOI']
            if not doi:
                raise ValueError
    #        print('doi')
        except (KeyError, ValueError):
            try:
                doi = summaryRecord['ArticleIds']['doi']
                if not doi:
                    raise ValueError
    #            print('article id doi')
            except (KeyError, ValueError):
                try:
                    doi = summaryRecord['ELocationID'].strip('doi: ')
                    if not doi:
                        raise ValueError
    #                print('elocation id doi')
                except (KeyError, ValueError):
    #                print('noDoi')
                    pass
        try:
    #        print(len(doi))
            if doi:
                return doi
            else:
                raise ValueError
        except NameError:
            raise ValueError



    def pmidTodoi (self, pmidsListOrCSVfile, doisCSVfile):
        # reading the pmids
    #    pmidsListOrCSVfile = 'related scores'
    #    doisCSVfile = 'dois'
        pmids = pmidsListOrCSVfile
        if type(pmids) is str:
            pmids = pd.read_csv(pmids+'.csv')
            pmids = set(pmids['PMID'])
            pmids = {str(pmid) for pmid in pmids}

        # retreive doi for each pmid
        try:
            unicodeErrorpmids = set(self.csvTolist('unicode Error pmids'))
            notFoundpmids = set(self.csvTolist('not Found pmids'))
            notFounddois = set(self.csvTolist('not Found dois'))
            pmid_doi = self.csvTodict(doisCSVfile)
            dois = set(pmid_doi.keys())
        except FileNotFoundError:
            unicodeErrorpmids = set()
            notFoundpmids = set()
            notFounddois = set()
            pmid_doi = dict()
            dois = set()

        pmids = pmids - unicodeErrorpmids - notFoundpmids - notFounddois - dois
        unicodeErrorpmids = set()
        notFoundpmids = set()
        notFounddois = set()
        pmid_doi = dict()
        dois = set()
        print(len(pmids))
        print('fetching..')
        counter = 0
        for pmid in pmids:
            counter+=1
    #        print(pmid, counter)

            try:
                summaryHandle = Entrez.esummary(db='pubmed', id=pmid)
            except URLError:
                print('server time out')
                time.sleep(60)
                print('woke up from server time out')
            try:
                summaryRecord = Entrez.read(summaryHandle)
                try:
                    doi = self.pmidEntrezSummaryRecordTodoi(summaryRecord)
                    pmid_doi[pmid] = self.extractDOI(doi).lower()
                except ValueError:
                    notFounddois.add(pmid)
    #                print('no doi', pmid, counter)
            except (UnicodeDecodeError, UnicodeEncodeError):
                unicodeErrorpmids.add(pmid)
                print('unicode', pmid, counter)
            except RuntimeError as e:
                notFoundpmid = re.findall('\d+', str(e))[0]
                assert pmid == notFoundpmid
                notFoundpmids.add(pmid)
                print('no record', pmid, counter)

            if counter % 1000 == 0 or counter == len(pmids):
                self.dictTocsv(pmid_doi, doisCSVfile, ['PMID', 'DOI'])
                self.listTocsvRows(unicodeErrorpmids,'unicode Error pmids')
                self.listTocsvRows(notFoundpmids, 'not Found pmids')
                self.listTocsvRows(notFounddois, 'not Found dois')
                time.sleep(15)
                pmid_doi=dict()
                unicodeErrorpmids = set()
                notFoundpmids = set()
                notFounddois = set()
                print('woke up')                    


    def writeXMLWithHeaderFrom(self, fileNameOfHeaderFile, fileNameOfOutFile, tree=None):
        if tree is None:
            tree = ET.parse(fileNameOfHeaderFile)
        root = tree.getroot()
        with open(fileNameOfHeaderFile, 'r') as inf, open(fileNameOfOutFile, 'wb') as outf:
            for line in inf:
                if line.strip()== '<'+ root.tag +'>':
                    break
                outf.write(line.encode('utf-8'))
            tree.write(outf)

    def getExtension(self, filePathOrURL):
        ext = re.findall('.+\.([a-zA-Z]{3,10})"?$', filePathOrURL)
        if ext:
            return ext[-1]
        else:
            return None

    def decideExtension(self, dictionaryOfcontentTypesDispositionURL):
        try:
            contentType = dictionaryOfcontentTypesDispositionURL['contentType']
        except KeyError:
            contentType = None
        try:
            crContentType = dictionaryOfcontentTypesDispositionURL['crContentType']
        except KeyError:
            crContentType = None
        try:
            contentDisposition = dictionaryOfcontentTypesDispositionURL['contentDisposition']
        except KeyError:
            contentDisposition = None
        try:
            url = dictionaryOfcontentTypesDispositionURL['url']
        except KeyError:
            url = None

        if contentType:
            ext = contentType.split('/')[-1].split(';')[0]
    #        print('type')
        elif crContentType and crContentType != 'unspecified':
            ext = crContentType.split('/')[-1].split(';')[0]
    #        print('cr')
        elif contentDisposition:
            ext = self.getExtension(contentDisposition)
    #        print('dispo')
        elif self.getExtension(url):
            ext = self.getExtension(url)
    #        print('url')
        else:
            ext = None
    #        print('non')

        if 'api.elsevier' in url and (ext == 'plain' or ext == 'htm'or ext == 'html'):
            ext = 'txt'
        elif ext == 'plain' or ext == 'htm':
            ext = 'html'
        elif ext == 'msword':
            ext = 'doc'
        elif ext == 'octet-stream':
            ext = 'pdf'

        return ext


    def multireplace(self, string, replacements):
        """
        From: https://gist.githubusercontent.com/bgusach/a967e0587d6e01e889fd1d776c5f3729/raw/a2ac838179d453a1b9a028cfee9a66ff0e0157dc/multireplace.py
        Given a string and a replacement map, it returns the replaced string.

        :param str string: string to execute replacements on
        :param dict replacements: replacement dictionary {value to find: value to replace}
        :rtype: str

        """
        # Place longer ones first to keep shorter substrings from matching where the longer ones should take place
        # For instance given the replacements {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
        # 'hey ABC' and not 'hey ABc'
        substrs = sorted(replacements, key=len, reverse=True)

        # Create a big OR regex that matches any of the substrings to replace
        regexp = re.compile('|'.join(map(re.escape, substrs)))

        # For each match, look up the new string in the replacements
        return regexp.sub(lambda match: replacements[match.group(0)], string)

    def extractDOI(self, doiLinkOrTag):
        return doiLinkOrTag.split('org/')[-1]

    def createFile(self, doi, ext=None):
        notAllowedCharReplacements = { 
            '/':'-',
            '\\':'-',
            '>':'}',
            '<':'{',
            ':':'_',
            '?':'!',
            '"':'\'',
            '*':'+',
            '|':'$'
        }
        if ext:
            fileName = 'full text/doi3/'+ext+'/'+self.multireplace(doi, notAllowedCharReplacements)+ '.'+ext
        else:
            fileName = 'full text/doi3/other/'+self.multireplace(doi, notAllowedCharReplacements)
        directory = os.path.dirname(fileName)
        if not os.path.exists(directory):
            os.makedirs(directory)
        return fileName


    def getDomain(self, url):
        domain = re.findall('//(.+?)/', url)
        if domain:
            return domain[-1]
        else:
            return None


    def fixSpringerLinkURL(self, url):
        doi = re.findall('pdf/(.+)', url)[-1]
        doi = doi.replace('/', '%2F')
        newURL = 'https://link.springer.com/content/pdf/' + doi
        return newURL


    def getCambridgeURL(self, doi):
        doi = re.findall('10.1017/s(.+)', doi)[-1]
        newURL = 'https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S' + doi
        return newURL


    def download(self, doi, url, requestHeaders):
        downloadableViaRequests = False
        try:
            r = requests.head(url, headers=requestHeaders, allow_redirects=True)
            try:
                remainingHits = int(r.headers['CR-TDM-Rate-Limit-Remaining'])
                resetTime = int(r.headers['CR-TDM-Rate-Limit-Reset'][: -3])
                if remainingHits < 2:
                    sleepingTime = abs(resetTime - time.time())
    #                print('sleeping')
                    time.sleep(sleepingTime+1)
    #                print('woke up')
            except KeyError:
                pass
            if r.status_code == requests.codes.ok:
                downloadableViaRequests = True
        except (requests.exceptions.ConnectionError,
                requests.exceptions.TooManyRedirects) as e:
    #        print(doi, url)
    #        print(str(e))
            pass
        if downloadableViaRequests:
            self.downloadViaRequests(doi, url, requestHeaders)
    #        print('downloading via requests')
        else:
            self.downloadViaUrlLib(doi, url, requestHeaders)
    #        print('downloading via urllib')


    def downloadViaRequests(self, doi, url, requestHeaders):
        try:
            content = requests.get(url, headers=requestHeaders, allow_redirects=True) #, stream=True)
            try:
                contentType = content.headers['content-type']
            except KeyError:
                contentType = None
            try:
                contentDisposition = content.headers['content-disposition']
            except KeyError:
                contentDisposition = None

            d = {'contentType':contentType,
                 'contentDisposition':contentDisposition,
                 'url':url
            }

            content = content.content
            contentSize = int(sys.getsizeof(content))
            fileSize = 0
            ext = self.decideExtension(d)
            fileName = self.createFile(doi, ext)
        #    print(fileName)
            if os.path.exists(fileName):
                fileSize = os.path.getsize(fileName)
        #        print(contentSize, fileSize)
            if contentSize-33 > fileSize:
                with open(fileName, 'wb') as f:
        #            print('writingFile')
        #            for chunk in content.iter_content(chunk_size=1024):
        #                if chunk:
                    f.write(content)

            fileName = 'full text/doi/downloadable dois'
            directory = os.path.dirname(fileName)
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.writecsvRow('full text/doi/downloadable dois', [doi, url, ext], ['DOI', 'URL', 'EXT'])
        except (requests.exceptions.HTTPError, requests.exceptions.Timeout, requests.exceptions.SSLError, requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError) as e:
            fileName = 'full text/doi/bad urls'
            directory = os.path.dirname(fileName)
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.writecsvRow('full text/doi/bad urls', [doi, url, str(e)], ['DOI', 'URL', 'ERROR'])
        except:
            print(doi, url)
            raise


    def downloadViaUrlLib(self, doi, url, requestHeaders):
        try:
            req = Request(url)
            if requestHeaders is not None:
                k, v = next(iter(requestHeaders.items()))
                req.add_header(k, v)
            urlibContent = urlopen(req)
            try:
                contentType = urlibContent.info()['Content-Type']
            except KeyError:
                contentType = None

            d = {'contentType':contentType,
                 'url':url
            }
            content = urlibContent.read()
            contentSize = sys.getsizeof(content)
            fileSize = 0
            ext = self.decideExtension(d)
            fileName = self.createFile(doi, ext)
    #        print(fileName)
            if os.path.exists(fileName):
                fileSize = os.path.getsize(fileName)
    #            print(contentSize, fileSize)
            if contentSize-33 > fileSize:
                with open(fileName, 'wb') as f:
    #                print('writingFile')
    #                for chunk in content.iter_content(chunk_size=1024):
    #                    if chunk:
                    f.write(content)

            fileName = 'full text/doi/downloadable dois'
            directory = os.path.dirname(fileName)
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.writecsvRow('full text/doi/downloadable dois', [doi, url, ext], ['DOI', 'URL', 'EXT'])
    #                print('created bad file')
        except (HTTPError, TimeoutError, URLError, ConnectionResetError, CertificateError) as e:
            fileName = 'full text/doi/downloadable dois'
            directory = os.path.dirname(fileName)
            if not os.path.exists(directory):
                os.makedirs(directory)
            self.writecsvRow('full text/doi/bad urls', [doi, url, str(e)], ['DOI', 'URL', 'ERROR'])
        except:
            print(doi, url)
            raise


    def csvTolist(self, csvFileName):
        cells = list()
        with open(csvFileName+'.csv', 'r') as inf:
            r = csv.reader(inf)
            for row in r:
                row = list(filter(None,row))
                if len(row) == 1:
                    row = row[0]
                if row:
                    cells.append(row)
            if len(cells) == 1:
                cells = cells[0]
        return cells


    def csvTodict(self, csvFileName):
        mydict = defaultdict(set)
        with open(csvFileName+'.csv', 'r') as f:
            csvin = csv.reader(f)
            for row in csvin:
                k, v = row[0], row[1]
                mydict = self.addToDict(mydict, k, v)
        return mydict

    def checkHeaders(self, csvFileName, headers):
        if headers is None and not os.path.isfile(csvFileName+'.csv'):
            raise FileExistsError('We need the headers to crearte the file for the first time.')
        elif headers and not os.path.isfile(csvFileName+'.csv'):
            with open(csvFileName+'.csv', 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)

    def dictTocsv(self, mydict, csvFileName, headers=None):
        self.checkHeaders(csvFileName, headers)
        unicodeDict=dict()
        with open(csvFileName+'.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            for key, value in mydict.items():
                try:
                    if type(value) is set or type(value) is list:
                        for x in value:
                            writer.writerow([key, x])
                    else:
                        writer.writerow([key, value])
                except (UnicodeEncodeError, UnicodeDecodeError):
                    unicodeDict[key] = value
        with open(csvFileName+'.csv', 'ab') as f:
            writer = unicodecsv.writer(f, encoding='utf-8')
            for key, value in unicodeDict.items():
                if type(value) is set or type(value) is list:
                    for x in value:
                        writer.writerow([key, x])
                else:
                    writer.writerow([key, value])

    def listTocsvCols(self, li, fileName):
        containsUnicode=False
        with open(fileName+'.csv', 'a', newline='') as csvFile:
            liFile = csv.writer(csvFile, quoting=csv.QUOTE_ALL)
            try:
                liFile.writerow(li)
            except (UnicodeEncodeError, UnicodeDecodeError):
                containsUnicode=True
        if containsUnicode:
            with open(fileName+'.csv', 'ab') as csvFile:
                liFile = unicodecsv.writer(csvFile, encoding='utf-8', quoting=csv.QUOTE_ALL)
                liFile.writerow(li)


    def listTocsvRows(self, li, fileName):
        containsUnicode=False
        directory = os.path.dirname(fileName)
        if '/' in directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(fileName+'.csv', 'a', newline='') as csvFile:
            liFile = csv.writer(csvFile, quoting=csv.QUOTE_ALL)
            try:
                liFile.writerows([i] for i in li)
            except (UnicodeEncodeError, UnicodeDecodeError):
                containsUnicode=True
        if containsUnicode:
            with open(fileName+'.csv', 'ab') as csvFile:
                liFile = unicodecsv.writer(csvFile, encoding='utf-8', quoting=csv.QUOTE_ALL)
                liFile.writerows([i] for i in li)


    def writecsvRow(self, fileName, row, headers=None):
        self.checkHeaders(fileName, headers)
        containsUnicode=False
        with open(fileName + '.csv', 'a', newline='') as f:
            csvout = csv.writer(f)
            try:
                csvout.writerow(row)
            except (UnicodeEncodeError, UnicodeDecodeError):
                containsUnicode=True
        if containsUnicode:
            with open(fileName+'.csv', 'ab') as f:
              csvout = unicodecsv.writer(f, encoding='utf-8')
              csvout.writerow(row)


    def addToDict(self, mydict, k, v):
        d = defaultdict(set, mydict)
        d[k] |= {v}
        return d


    def removeKey(self, d, keys):
        r = dict(d)
        try:
            for key in keys:
                try:
                    del r[key]
                except KeyError:
                    pass
        except TypeError:
            del r[keys]
        except KeyError:
            pass
        return r
