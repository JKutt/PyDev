import PyPDF2 

filename = "C:\\Users\\johnk\\devProjects\\Python\\PyDev\\pdf\\Nov 26 2018_Daily Activity Report.pdf"
#open allows you to read the file
pdfFileObj = open(filename,'rb')
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
splitted = sub_text.split(" ")
print(splitted)