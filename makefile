test:
	gcc -o test source/test.c -lm -lncurses -g
gravitas: 
	g++ -o binaries/gravitas-linux-x86_64 source/gravitas.cpp -lGL -lGLEW -lGLU -lGL -lglfw -lm -lXrandr -lXi -lXxf86vm -lpthread -g
clean: 
	rm binaries/*
