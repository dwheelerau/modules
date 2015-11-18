#!/usr/bin/python
import smtplib  

class Send_email:
    def __init__(self,toaddrs,msg):
        '''sends an email, username is just string before @gmail *req's gmail'''
        self.fromaddr = "dwheelerau.work@gmail.com"
        self.toaddrs = toaddrs
        self.pwd = "Pubcrl75#!"
        self.msf = msg
        self.user_name = "dwheelerau.work"

    def send_email(self):
        server = smtplib.SMTP("smtp.gmail.com:587")
        server.ehlo()
        server.starttls()
        server.ehlo()
        try:
            server.login(self.user_name,self.pwd) 
            server.sendmail(self.fromaddr, self.toaddrs, self.msf)  
        except smtplib.SMTPAuthenticationError:
            print "email-not sent - login_fail!"
        server.quit()  
     
