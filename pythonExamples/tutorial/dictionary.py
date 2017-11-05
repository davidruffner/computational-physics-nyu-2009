#dictionary
phonebook={'Andrew':8,\
           'Emily Everett':6,'Pete':7,\
           'Lewis':1}
print phonebook
print "here are the numbers", phonebook.values()

phonebook['G']=12
print "added a number, now"
print phonebook.values()

del phonebook['Lewis']

print "deleted an entry"
print phonebook.values()

if phonebook.has_key('Andrew'):
    print "Andrew is in the dict. She is", \
          phonebook['Andrew']
else:
    print "Andrew is not in the dict"

print "They are in the phonebook:"
print phonebook.keys()

keys = phonebook.keys()

keys.sort()
print keys

values = phonebook.values()
print values
values.sort()
print values

print "The dictionary has", \
      len(phonebook), "entries in it"
