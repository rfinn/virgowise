#PURPOSE:
#To recreate errors produced by galfit so that they can be addressed in fitting of WISE imaging
#GOAL:
#Define the errors so that g.py can recognize when they happen
import os


#####################################Stuff for 'crashing', when galfit blows up#######################
#Do this if corrupting a WISE file doesn't make self.error = 1

###########

#From what I can understand from the manual, self.error will be 1 when galfit crashes, the way we're running it (page 26 of https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf, top, right, 3rd paragraph down)

#I added a break in unwise_rungalfit_dev in process_list  to provide the function of reacting if this error is a 1, which works. Currently beak just stops process_list fromcontinuing to loop and exits.

#The goal will be to (I think) note when this happens but proceeed, so the galaxies on which galfit does this can be re-evaluated. This should be pretty simple, we could create a file to which we append the galaxy names which error = 1 onto

###########


#probably don't need this next part
#if self.error = 1:
    #Highlight for re-evaluation. Caused by galfit not being able to use inputs 

#errtest = os.system('galfit'+self.galfile >> 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'error_test.txt') #grabs galfit output
#contents = #make string of 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'error_test.txt' 
#find_bomb = '@@@'
#contents.find(find_bomb)


#####################################Stuff for a failure to converge##################################  
#check = self.xc
#if check.find(*)=-1
    #Then, don't use this fit. Highlight for re-evaluation
#It seems like this is also noted and stored in self.error. When running unwise_galfit_dev right now it always breaks when a warning for possible numerical problems with output parameters happens. 
    
    
