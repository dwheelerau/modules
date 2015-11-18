#!/usr/bin/python
'''a script for doing the stats'''
from statlib import stats ##requires statlib-1.1.tar.gz
from math import log, exp

def mann_white(venom_vals_list,others_gp_vals):
    '''compares one mean to another are they signif different
    takes a two list of floats for each group
    Mann-... U test for significant difference between
    two groups\
    >>> 
    >>>mann_white(list_venom_dnds_as_floats,list_comparable_gp_vals)
    ...Mean venom dnds = 0.470993103448, mean others 0.292471404883.
    ...Mann-whitney test U = 205437.0 with a p(one tailed) = 2.42740449563e-06
    >>> 
    One-sided p-value assuming a asymptotic normal distribution.
    Notes

    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is significant
    if the u-obtained is LESS THAN or equal to the critical value of U.

    This test corrects for ties and by default uses a continuity correction
    . The reported p-value is for a one-sided hypothesis, to get the two-sided
    p-value multiply the returned p-value by 2.'''
    test_gp = [vals for vals in venom_vals_list] #redundent but helps understanding
    other_gp = [vals for vals in others_gp_vals]
    result = stats.mannwhitneyu(test_gp,other_gp)

    test_mean = stats.mean(test_gp)
    other_mean = stats.mean(other_gp)
    return "Mean venom dnds = %s, mean others %s.Mann-whitney test U = %s with a p(one tailed) = %s"\
      %(test_mean,other_mean,result[0],result[1])

def mann_white_alt(venom_vals_list,others_gp_vals):
    '''compares one mean to another are they signif different
    takes a two list of floats for each group
    Mann-... U test for significant difference between
    two groups\
    >>> 
    >>>mann_white(list_venom_dnds_as_floats,list_comparable_gp_vals)
    ...Mean venom dnds = 0.470993103448, mean others 0.292471404883.
    ...Mann-whitney test U = 205437.0 with a p(one tailed) = 2.42740449563e-06
    >>> 
    One-sided p-value assuming a asymptotic normal distribution.
    Notes

    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is significant
    if the u-obtained is LESS THAN or equal to the critical value of U.

    This test corrects for ties and by default uses a continuity correction
    . The reported p-value is for a one-sided hypothesis, to get the two-sided
    p-value multiply the returned p-value by 2.'''
    test_gp = [vals for vals in venom_vals_list] #redundent but helps understanding
    other_gp = [vals for vals in others_gp_vals]
    result = stats.mannwhitneyu(test_gp,other_gp)

    test_mean = stats.mean(test_gp)
    other_mean = stats.mean(other_gp)
    return "Mean venom dnds = %s, mean others %s.\nMann-whitney test U = %s with a p(one tailed) = %s"\
      %(test_mean,other_mean,result[0],result[1])
def raw_bincoeff(n, r):
    '''Calculates the binomial coefficient'''

    if r < n - r:
        r = n - r
    x = 1
    for i in range(n, r, -1):
        x *= i
    for i in range(n - r, 1, -1):
        x /= i
    return x

#MODIFIED FROM DHYPER.PY REEF SOFTARE
#The hypergeometric distribution has the form,
#
#   prob(k) =  choose(n1,t) choose(n2, t-k) / choose(n1+n2,t)
#
#where choose(a,b) = a!/(b!(a-b)!) 
#
#n1 + n2 is the total population (tagged plus untagged)
#n1 is the tagged population
#t is the number of samples taken (without replacement)
#k is the number of tagged samples found
def hypergeometric(k, n1, n2, t):
    '''n1=succ pop tot,k num succ sel, n2 num fail, t sampels taken
    calculates probably of selecting k things from a group(n1+n2) with
    n1 things available and n2 things in the other gp ie.
    50 marbles, 10 red, 40 black, what change get 5 red from 10 tries..
    k=5,n1=10,n2=40,t=10 - used when prior selectiong effects future
    if independent ie coin toss then would used binomial dist (wit replacment)
    see www.stattrek.com for some good advice here or wiki
    '''
    #If the value of t is wrong
    if t > n1 + n2:
        t = n1 + n2
    if k > n1 or k > t:
        return 0
    elif t > n2 and ((k + n2) < t):
        return 0
    else:
        c1 = log(raw_bincoeff(n1,k))
        c2 = log(raw_bincoeff(n2, t - k))
        c3 = log(raw_bincoeff(n1 + n2 ,t))
    
        return exp(c1 + c2 - c3)
