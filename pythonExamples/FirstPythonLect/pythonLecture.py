import sys

list = [1,2,3,11,3,2,111,9]
print list

print sorted(list)
list.sort()

print list.sort()
print list

for x in reversed(sorted(list)):
    print x

filename = "testFile.txt"

try:
    testFile = open(filename)
except IOError, e:#the error message is saved in e
    print "Cant open the file, here is the error message"
    print e
    sys.exit()
    
for line in testFile:
    if "#" in line:
        continue
    print line

print "done"
