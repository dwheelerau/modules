import re
from decimal import *
import urllib

class Share:
    '''a class for a share ticker and number and access to price by web'''
    def __init__(self,ticker,number,purchase_price):
        self.ticker = ticker
        self.number = int(number)
        self.purchase_price = Decimal(purchase_price)
        self.value = self.get_value()
        
    def get_value(self):#,ticker):
        base_url = 'http://finance.google.com/finance?q='
        content = urllib.urlopen(base_url + self.ticker).read()
        m = re.search('span id=".*?>(.*?)<', content)
        if m:
            quote = m.group(1)
        else:# m==None:
            print 'no quote available for: ' ,self.ticker
            return "0.0" #needs to return a string digit for decimal function
        return Decimal(quote)
            
