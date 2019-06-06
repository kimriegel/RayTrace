import time
import numpy
## List
#def square(nums):
#    result = []
#    for i in nums:
#        result.append(i*i)
#    return result
#
#nums = square(range(1,6))
#print(nums)

## Generator
#
#def square(nums):
#    for i in nums:
#        yield(i*i)      # What makes it a generator
#
#nums = square(range(1,6))
#print(nums)             #prints the object in memory
#print(next(nums))       #Prints first point in list
#
#for n in nums:
#    print(n) 
#                        # Will run 4 times if accessed outside of this loop, 5 if not
#                        #  It will iterate through everything once and only once

## List comprehension   
#Cleaner version of part 1
#nums = [i*i for i in range(1,6)]
#for n in nums:
#    print(n) 

## Generator with list comprehension
## Generator, but cleaner
#nums = (i*i for i in range(1,6))    #Changing from brackets to parenthesis makes this work differently
#print (nums)
#for n in nums:
#    print(n) 

    # Checking time

#def square_list(nums):
#    result = []
#    for i in nums:
#        result.append(i*i)
#    return result
#
#
#def square_gen(nums):
#    for i in nums:
#        yield(i*i)      # What makes it a generator
#
#t= time.time()
#nums = square_list(range(1,100001))
#print(time.time()-t)
#
#t= time.time()
#nums = square_gen(range(1,100001))
#print(time.time()-t)

#print(nums)             #prints the object in memory
#print(next(nums))       #Prints first point in list
#
#for n in nums:
#    print(n) 
#                        # Will run 4 times if accessed outside of this loop, 5 if not
#                        #  It will iterate through everything once and only once

# Generator

def carpet2(rayrange):
    for i in rayrange:
        yield (i,i*i)      # What makes it a generator

rays = carpet2(range(6))
print(rays)             #prints the object in memory
#print(next(nums))       #Prints first point in list

for veci in rays:
    print(veci) 
                        # Will run 4 times if accessed outside of this loop, 5 if not
                        #  It will iterate through everything once and only once
