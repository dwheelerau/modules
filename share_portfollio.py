#from decimal import *
from Share import *

class Portfolio:
    '''a class of share holdings takes a dictonary of share:number'''
    def __init__(self,share_protfolio_dic):
        self.share_portfolio_dic = share_protfolio_dic
        self.share_list = []
        self.old_value = []
        for share in self.share_portfolio_dic:
            ticker,number,purchase_price,old_value = share[0],share[1],\
                                                     share[2],share[3]
            self.old_value.append(share[3])
            self.share_list.append(Share(ticker,number,purchase_price))#,old_value)) #send ticker to class
        self.old_value = sum([Decimal(a) for a in self.old_value])
        #print sum(self.old_value)
        self.holding_value = self.get_holding_value()

    def get_holding_value(self):
        bank = Decimal(0)
        for share in self.share_list:
            price = Decimal(str(share.value))
            holding=int(share.number)
            tot_value = price*holding
            bank = bank+tot_value
        the_difference = str(100*((bank-self.old_value)/self.old_value))
        rounded_diff = the_difference[:the_difference.find(".")+2]
        return bank, rounded_diff
    
    def add_to_holding(self, ticker,number,price_per_share):#add purchase price
        #add new shares or increase decerease
        share_found = "No"
        for n in self.share_list:
            if n.ticker == ticker.upper():
                cost = (n.purchase_price*n.number) + (Decimal(price_per_share)*number)
                n.number=n.number + number
                n.purchase_price = cost/n.number
                share_found = "yes"
        if share_found=="No":
            #old_value = number*price_per_share
            #self.holding_value=self.holding_value+old_value
            self.share_list.append(Share(ticker.upper(),number,price_per_share))#,old_value))
            
    def sell_shares(self,ticker,number="All"):
        if number=="All": #default sell all shares
            for n in self.share_list:
                if n.ticker==ticker.upper():
                    self.share_list.remove(n)
        else:
            for n in self.share_list:
                if n.ticker==ticker.upper():
                    n.number=n.number-number
                    
    def tweet_it(self):#could put this in alarm class? Alarm may have to inherit stuff?
        json=0
        usr,pwd = "dwheelerau",raw_input("please enter password for dwheelerau: ")
        constring = "http://%s:%s@twitter.com/statuses/update.json" % (usr,pwd)
        tweet = ""
        for n in self.share_list:
            tweet="%s %s "%(tweet,n.ticker)#get it to send value I think?
        tweet = urllib.urlencode({"status":tweet})
        f = urllib.urlopen(constring, tweet)
        if json ==1:
            print f.read()
        else:
            print "tweet sent"
            
