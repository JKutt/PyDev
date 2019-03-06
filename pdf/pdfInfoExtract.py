import PyPDF2 
from os.path import isfile, join, isdir
from os import listdir

# filename = ["C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\pdfs\\Nov 26 2018_Daily Activity Report.pdf"]
dir_stuff = "C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\pdfs\\"
dir_path = [f for f in listdir(dir_stuff) if isfile(join(dir_stuff, f))]     # get files
output_file = "C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\line_totals.txt"

for i in range(len(dir_path)):
	#open allows you to read the file
	pdfFileObj = open(dir_stuff + dir_path[i],'rb')
	#The pdfReader variable is a readable object that will be parsed
	pdfReader = PyPDF2.PdfFileReader(pdfFileObj)
	#discerning the number of pages will allow us to parse through all #the pages
	num_pages = pdfReader.numPages
	print(num_pages)
	pageObj = pdfReader.getPage(0)
	text = pageObj.extractText()
	# lets start with productivity
	start_productivity = text.find("1.7") + 3
	end_productivity = text.find("Environmental")
	sub_text = text[start_productivity:end_productivity]
	out_file = open(output_file, "w+")
	out_file.write('%s' % sub_text)
	print(sub_text)
# splitted = sub_text.split(" ")
# print(splitted)