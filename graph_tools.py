#!/usr/bin/python
import matplotlib.pyplot as plt

'''see matlabplot.txt for help'''
def draw_line_plot(x_list,y_list,y_err_list, outfile="",x_lab="",y_lab="",title=""):
    '''Draws sing/multi lines on a single graph, unpacks list of x points and
    list of y points and list of yerrors:
    usage(multi = ([[1,2,3],[1,2,3]],[[1,2,3],[3,4,5]],[[0.1,0.2,0.3],[0.1,0.2,0.3]])
    single = ([1,2,3],[1,2,3] etc
    output format controlled by by *.jpg, *png etc. empty output
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colours = ["r","g"]
    labels = ["red","green"]
    for n in range(len(x_list)):
        x = x_list[n]
        y = y_list[n]
        err = y_err_list[n]
        ax.plot(x,y)
        ax.errorbar(x,y, yerr=err, fmt='o')
        
        #add lines to indicate catagory 1,2,3 median vals of cpgoe for each shown
        #         ax.axhline(y=0.98, xmin=0, xmax=1,color="g",label="test")
        #ax.axhline(y=1.5, xmin=0, xmax=1,color="b")
        #ax.axhline(y=1.88, xmin=0, xmax=1,color="r")

        #add xlab
        ax.set_xlabel(x_lab)
        ax.set_ylabel(y_lab)
        #tit_str = "%s average CPGOE in a non-overlapping window size = %s genes"%(scaf,WINDOW)
        ax.set_title(title)
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()  

#def line_plot(x_list,y_list):
#    pass
#
#def line_plot_error(x_list,y_list,y_err_list, outfile="",x_lab="",y_lab="",title=""):
#    ''' draw graph, outfile default will print, output empty will print to screen
#    output format controled by *.jpg, *png etc''' 
#    # make an object
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(x_list,y_list)
#    ax.errorbar(x_list,y_list, yerr=y_err_list, fmt='o')
#    
#    #add lines to indicate catagory 1,2,3 median vals of cpgoe for each shown
#    #ax.axhline(y=0.98, xmin=0, xmax=1,color="g",label="test")
#    #ax.axhline(y=1.5, xmin=0, xmax=1,color="b")
#    #ax.axhline(y=1.88, xmin=0, xmax=1,color="r")
#
#    #add xlab
#    ax.set_xlabel(x_lab)
#    ax.set_ylabel(y_lab)
#    #tit_str = "%s average CPGOE in a non-overlapping window size = %s genes"%(scaf,WINDOW)
#    ax.set_title(title)
#    if outfile:
#        plt.savefig(outfile)
#    else:
#        plt.show()
    

