list1 = range(50)
list2 = list()
for x in list1:
    list2.append(x**2)

print list2

'''or more simply'''
list3 = [x**3 for x in list1]
print list3

'''more options'''
list4 = [x**3 for x in list1 if x % 2 ==0]
print list4
