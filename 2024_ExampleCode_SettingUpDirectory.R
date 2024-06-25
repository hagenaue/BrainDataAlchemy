#Sample code for setting up an R project directory
#Megan Hagenauer, May 8 2024

#################### 

#0. We'll start with a basic tour:

#Files
#Console
#Script editor (.R file)
##Code short cuts: run code, comment/uncomment

#Tools: Code packages

#Environment:
## Naming: objects need to start with letters, can contain numbers, _, ., but no spaces. Lowercase recommended.
## Objects created with code *are not saved as part of the script (.R) file*
## Session: Working Directory, Workspace (.RData file)

#Other: packages, help, plots

#Useful reference:
# https://r4ds.hadley.nz/workflow-scripts.html
#Also has a nice section on *file naming*


####################

#1. Tell R where you want your home directory to be (i.e., where your R project folders will be located in general - not just this project):

#The tilde is shorthand for your current home directory, as defined in R-Studio. You can see the default using this code: 
path.expand("~/")
#Mine is:
#[1] "/Users/hagenaue/"

#Example: I would like to change my home directory to my Documents Folder: 
#Note - if you are on a PC instead of a mac, you may have the direction of slashes in your file path reversed from / to \, but R doesn't understand \ so you have to either use / or \\
Sys.setenv(HOME = "/Users/hagenaue/Documents")

#If I wanted to revert to the original home directory, I would just set it back to the original path:
#Sys.setenv(HOME = "/Users/hagenaue")

#Double-checking that my home directory is now the Documents folder:
path.expand("~/")
#[1] ""/Users/hagenaue/Documents/""

#Double-checking (again) that my home directory is my Documents folder by looking at whether it has the same contents as my Documents folder:
list.files("~/")

###########################

#2. Create a folder for all files (input, code, output, coding environment) for your project on your computer and then set the working directory to that folder.

#Creating a directory (folder) for your project files:
dir.create("~/Summer2024_MetaAnalysisProject")

#Now I'm going to tell R that I want to use my new folder as the working directory for my project:
#A working directory is where R looks for files that you ask it to load, and where it will put any files that you ask it to save
setwd("~/Summer2024_MetaAnalysisProject")

#You can also do this in Rstudio via the GUI drop down menu (Session-> Set Working Directory-> Choose Directory)
#but if you do, you should save the code that appears in the Console so that you can easily navigate back.

#If you want to double-check which working directory you are in:
getwd()


#Let's double-check that there isn't anything currently in our working directory:
list.files()
#character(0)
#Nope :)

#And make some folders for all of our project components:

dir.create(path="~/Summer2024_MetaAnalysisProject/R_Input")
dir.create(path="~/Summer2024_MetaAnalysisProject/R_Output")
dir.create(path="~/Summer2024_MetaAnalysisProject/R_Figures")

#Now we can see those folders:
list.files()
#[1] "R_Figures" "R_Input"   "R_Output"


#################### 

#3. Creating a reproducible coding environment: Automatically tracking the code packages (and their versions!) that you are using for your project:

#Installing a code package that let's you track and reproduce the versions of the code packages that you will be using for your analysis:
install.packages("packrat")  

#A packrat project has its own private package library. 
#Any packages you install from inside a packrat project are only available to that project.
#packages you install outside of the project are not available to the project.

#Here's an overview of how packrat works:
# https://rstudio.github.io/packrat/walkthrough.html

# To track the versions of the packages that you are using for your project, link to the project directory (copied above) to use Packrat with packrat::init:
packrat::init()

#It will then restart the R session

