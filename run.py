import os 
import time

is_video = False
# is_video = True



if(is_video):
    # os.system("make clean")
    os.system("make")
    # os.system("clear")
    t1 = time.time()
    os.system("./main")
    t2 = time.time()
    print("Time: " + str(t2 - t1) + "s")
    os.system("python Pylib/genGIF.py")
    # os.system("open output.tga")
else :
    # os.system("make clean")
    os.system("make")
    # os.system("clear")
    t1 = time.time()
    os.system("./main")
    t2 = time.time()
    print("Time: " + str(t2 - t1) + "s")
    # os.system("python genGIF.py")
    os.system("open output.tga")



