import matplotlib.pyplot as plt
import numpy as np
import Jtools as JT

# plot second signal
##fileName = "C:\\Users\\johnk\\Projects\\ts1.xyz"
##fileName2 = "C:\\Users\\johnk\\Projects\\ts2.xyz"
##fileName3 = "C:\\Users\\johnk\\Projects\\ts5.xyz"
fileName = "C:\\Users\\johnk\\Projects\\test-02.xyz"
fileName2 = "C:\\Users\\johnk\\Projects\\test-03.xyz"
# fileName3 = "C:\\Users\\johnk\\Projects\\synctest0.xyz"

lines = 0
splarray = " "

# open the file
text_file = open(fileName,"r")

# determine how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()
# intiate numpy variables
X1 = np.zeros((lines))
Z1 = np.zeros((lines))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName,"r")
# retrieve and assign the data
for line in text_file:
        spl = line.split()
        X1[cnt] = float(spl[0])
        Z1[cnt] = float(spl[1])
        cnt += 1

text_file.close()

# open the second file
text_file = open(fileName2,"r")
lines =0
# determine how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()
# intiate numpy variables
X2 = np.zeros((lines))
Z2 = np.zeros((lines))

# initiate data count
cnt = 0
# initiate header count
cnthdr = 0
# open file again
text_file = open(fileName2,"r")
# retrieve and assign the data
for line in text_file:
        spl = line.split()
        X2[cnt] = float(spl[0])
        Z2[cnt] = float(spl[1])
        cnt += 1

text_file.close()
# print(X1[0])
# X1 = X1 + ((1. / 150. * 4.) * (1. / 60. * 1. / 60. * 1. / 24.))
# print(X1[0])
# open the third file
##text_file = open(fileName3,"r")
##lines =0;
###determin how many lines in the file
##while text_file.readline():
##        lines += 1
##print(lines)
##text_file.close()
### intiate numpy variables
##X3 = np.zeros((lines))
##Z3 = np.zeros((lines))
##
### initiate data count
##cnt = 0
### initiate header count
##cnthdr = 0
### open file again
##text_file = open(fileName3,"r")
### retrieve and assign the data
##for line in text_file:
##        spl = line.split()
##        X3[cnt] = float(spl[0])
##        Z3[cnt] = float(spl[1])
##        cnt += 1;

##X2 = X2+(160e-9)
##text_file.close()
##print(X2.size)
##sep = np.zeros(X1.size-1)
##cnt77 = 0
##cnt78 = 0
##for r in range(X1.size - 1):
##        sep[r] = (X1[r+1]*1e9)-(X1[r]*1e9)
##        if sep[r] == 77:
##                cnt77+=1
##        if sep[r] == 78:
##                cnt78 +=1
        
##print((sep))
##print("mean of timp sep")
##print((cnt77),'%0.2f')
##print((cnt78),'%0.2f')
##print(np.max(sep),np.min(sep))
##
##print("now lets check samples")
##print("Time Span: ",(np.floor((X2[X2.size-1] - X2[0])*86400)))
##print("numSamp Theoretical: ",(np.floor((X2[X2.size-1] - X2[0])*86400)*150.0))
##print("actual: ",X2.size)
##
####ts1trimX1 = X1[marpos-1:(X1.size-1)]
####ts1trimZ1 = Z1[marpos-1:(X1.size-1)]
####ts2trimX2 = X2[0:marend]
####ts2trimZ2 = Z2[0:marend]
####        print(i)
##
##marpos = 0
##marend = 0
### lets look for starting points
##for r in range(X1.size):
##        if X1[r] >= X2[0]:
##                marpos = r
##                break
##for r in range(X2.size):
##        if X2[r] >= X1[X1.size-1]:
##                marend = r
##                break
##
##print(marpos,marend)
##print(X1[marpos-1],'%0.9f')
##print(X2[0],'%0.9f')
##print(X1[marpos],'%0.9f')

##test = np.linspace(0, 239, num=35850)
##test = (test/86400.0)+X1[0]
##val = np.zeros((test.size))
##print(test.size)
##for i in range(test.size):
##        val[i] = JT.interpts(X1,Z1,test[i])

##fig = plt.figure()
##ax = fig.add_subplot(111)
##plt.plot(X2[0:100],Z2[0:100],'r.-')
##for xy in zip(X2[0:100], Z2[0:100]):                                       # <--
##    ax.annotate('(%0.9f, %s)' % xy, xy=xy, textcoords='data') # <--
##    plt.plot(X2[0:100],Z2[0:100])
##
##plt.plot(X1[300:420],Z1[300:420]-175,'.-')
##for xy in zip(X1[300:420], Z1[300:420]-175):                                       # <--
##    ax.annotate('(%0.9f, %s)' % xy, xy=xy, textcoords='data') # <--
##    plt.plot(X1[300:420],Z1[300:420]-175)
##plt.grid()
####plt.show()

plt.plot((X2), Z2, 'o-k')
plt.plot((X1), Z1, 'o-r')
plt.grid()
##plt.plot((ts2trimX2+160e-9),ts2trimZ2,'k.-')
##plt.plot(ts2trimZ2,'k')
##plt.plot(test,val,'r')
####plt.ylabel('signal')
######plt.hold(True)
#####plt.plot(X1[0:17099],Z2,'red')
##plt.plot(test,val,'g')
##plt.plot(X2,Z2,'r')
##plt.plot(ts1trimZ1,'r')

##plt.plot(ts2trimZ2-ts1trimZ1,'g')
#####plt.plot(X1[2399],Z1[2399],'g^')
#####plt.plot(X1[7351],Z1[7351],'g^')
#####plt.plot(X2,Z2,'red')
plt.show()

print("Plotting Complete")
