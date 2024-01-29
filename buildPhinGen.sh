g++ -std=c++11 -Wno-deprecated -Wno-unknown-pragmas -D_REENTRANT -DNDEBUG -DSC_MEMORY_ALIGNMENT=1 -DSC_LINUX -Wno-deprecated -fPIC -I ../supercollider-master/include/plugin_interface/ -I ../supercollider-master/include/common/ -I ../supercollider-master/include/server/ -g -c phingen3.cpp

g++ -shared -g -o phingen3.so phingen3.o

cp Segments.sc ~/.local/share/SuperCollider/Extensions/
cp FrAmpSegs.so ~/.local/share/SuperCollider/Extensions/
