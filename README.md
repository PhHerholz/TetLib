# TETLIB

Clone this repository. Change to the repo's directory

	cd tetlib

Clone ligibl into the same folder

	git clone https://github.com/libigl/libigl.git

Create and move to build folder

	mkdir build && cd build 

Compile everything using make. 

	cmake .. && make -j 8

Or use the right generator for your project files if you want to use an IDE i.e.
	
	cmake .. -G"Xcode" 

Before running, make sure the ./data folder is one level above the executable (this might be different for different build environments).
