#TEXT ADVENTURE GAME

#the menu function:
def menu(list, question):
    for entry in list:
        print 1 + list.index(entry),
        print ") " + entry

    return input(question) - 1

#the inspect function
def inspect(choice,location):
    if choice == location:
        print ""
        print "You found a key!"
        print ""
        return 1
    else:
        print ""
        print "Nothing of interest here."
        print ""
        return 0


#Give the computer some basic information about the room:
items = ["pot plant","painting","vase","lampshade","shoe","door"]

#The key is in the vase (or entry number 2 in the list above):
keylocation = 2

#You haven't found the key:
keyfound = 0

loop = 1

#Give some introductary text:
print "Last night you went to sleep in the comfort of your own home."

print "Now, you find yourself locked in a room. You don't know how"
print "you got there, or what time it is. In the room you can see"
print len(items), "things:"
for x in items:
    print x
print ""
print "The door is locked. Could there be a key somewhere?"
#Get your menu working, and the program running until you find the key:
while loop == 1:
    keyfound = inspect(menu(items,"What do you want to inspect? "),keylocation)
    if keyfound == 1:
        print "You put the key in the lock of the door, and turn it. It opens!"
        loop = 0


# Remember that a backslash continues
# the code on the next line

print "Light floods into the room as \
you open the door to your freedom."

