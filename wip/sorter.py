grid = [(93.4428213 ,28.8397178, 0.151),
        (68.5832    ,28.5998   , 0.151),
        (62.5832    ,28.5998   , 7.9423),
        (-2.40793   ,31.5003401, 0.151),
        (75.11005   ,28.4945787, 0.151)]

# Create short functions to parse list
# These are what parser will pay attention to when sorting
xvalue = lambda obj:obj[0]
yvalue = lambda obj:obj[1]
zvalue = lambda obj:obj[2]

# Sort list 
print(grid)

# A band is the strip of the environment to be checked, these are its values 
bando = 10 # Band length
bandmin = -10
bandmax = 0 + bando 
# We will be separating by x.
# X was chosen because it's easier to visually differentiate, but we could use any coordinate.
for i in range(10):
    band = [(x,y,z) for (x,y,z) in grid if ((bandmin<x)and(x<bandmax)) ]   # gives syntax error if for loop is excluded
    print(band)
    bandmin += bando
    bandmax += bando

# We can sort either the whole grid or the bands, both will give the same end result, but the single grid will be faster
grid.sort(key=xvalue)
bandmin = -10
bandmax = 0 + bando
for i in range(10):
    band = [(x,y,z) for (x,y,z) in grid if ((bandmin<x)and(x<bandmax)) ]   # gives syntax error if for loop is excluded
    #band.sort(key=xvalue)  # Alternatively this, but it's slower
    print(band)
    bandmin += bando
    bandmax += bando