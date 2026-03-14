all: main 
	echo "Hello World"

main: compute compute.cpp
	g++ compute.cpp -Lbuild -lcompute -I./src/ -I./include -o ./build/compute -lz

compute: 
	nvcc --shared -Xcompiler -fPIC -o ./build/libcompute.so ./src/compute.cu
	cp ./build/libcompute.so ./libcompute.so

transfer:
	g++ transfer.cpp -I./include -o ./build/transfer


