import PyPDF2 
from os.path import isfile, join, isdir
from os import listdir

# filename = ["C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\pdfs\\Nov 26 2018_Daily Activity Report.pdf"]
dir_stuff = "C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\pdfs\\"
dir_path = [f for f in listdir(dir_stuff) if isfile(join(dir_stuff, f))]     # get files
output_file = "C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\line_totals.txt"

out_file = open(output_file, "w+")
for i in range(len(dir_path)):
	#open allows you to read the file
	pdfFileObj = open(dir_stuff + dir_path[i],'rb')
	print("reading file {0}".format(dir_stuff + dir_path[i]))
	#The pdfReader variable is a readable object that will be parsed
	pdfReader = PyPDF2.PdfFileReader(pdfFileObj)
	#discerning the number of pages will allow us to parse through all #the pages
	num_pages = pdfReader.numPages
	pageObj = pdfReader.getPage(0)
	text = pageObj.extractText()

	# find date
	date_1_start = text.find("Current Date") + 15
	date_1_end = text.find("1.3") - 21
	out_file.write('%s -- ' % text[date_1_start:date_1_end])

	# lets start with productivity
	start_productivity = text.find("Productivity - Number of current injections and nominal current injection spacing") + 84
	end_productivity = text.find("Environmental")
	print(start_productivity , end_productivity)
	if (start_productivity == -1):
		pageObj = pdfReader.getPage(1)
		text = pageObj.extractText()
		start_productivity = text.find("Productivity - Number of current injections and nominal current injection spacing") + 84
		end_productivity = text.find("Environmental")
		sub_text = text[start_productivity:end_productivity]
		out_file.write('%s\n\n' % sub_text)
	elif (end_productivity == -1):
		if (start_productivity > 500):
			start_productivity = text.find("Productivity - Number of current injections and nominal current injection spacing") + 84
			sub_text = text[start_productivity:end_productivity]
			pageObj2 = pdfReader.getPage(1)
			text2 = pageObj2.extractText()
			end_productivity = text2.find("Environmental")
			sub_text_2 = text2[0:end_productivity]
			if (sub_text == "Page 1 of "):
				sub_text = sub_text_2
			else:
				sub_text = sub_text + sub_text_2
			out_file.write('%s\n\n' % sub_text)
		else:
			pageObj2 = pdfReader.getPage(1)
			text2 = pageObj2.extractText()
			start_productivity = text2.find("Productivity - Number of current injections and nominal current injection spacing") + 84
			end_productivity = text2.find("Environmental")
			sub_text = text2[start_productivity:end_productivity]
			out_file.write('%s\n\n' % sub_text)
	else:
		sub_text = text[start_productivity:end_productivity]
		out_file.write('%s\n\n' % sub_text)		
