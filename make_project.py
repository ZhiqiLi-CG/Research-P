import os
import sys
import platform

if platform.system().lower() == 'windows':
    sysid='win'
else:
    sysid='unix'

if __name__=='__main__':

    # first, update submodule
    cmd= "git.exe submodule update --progress --init --force"
    print(cmd)
    os.system(cmd)
    print("Finish Update Submodule")
    print("===============================================================")
    print("===============================================================")
    
    # build the visulization
    origin_root=os.path.join('.')
    build_root=os.path.join(origin_root,'build')
    if sysid=="win":
        cmd="cmake -S {} -B {} -G \"Visual Studio 17 2022\"".format(origin_root, build_root)
    else:
        cmd="cmake -S {} -B {} -DCMAKE_BUILD_TYPE=Release".format(origin_root, build_root)
    print(cmd)
    os.system(cmd)
    build_root=os.path.join(origin_root,'build')
    
    origin_root=os.path.join('.','Research-M')
    build_root=os.path.join(origin_root,'build')
    if sysid=="win":
        cmd="cmake -S {} -B {} -G \"Visual Studio 17 2022\"".format(origin_root, build_root)
    else:
        cmd="cmake -S {} -B {} -DCMAKE_BUILD_TYPE=Release".format(origin_root, build_root)
    print(cmd)
    os.system(cmd)