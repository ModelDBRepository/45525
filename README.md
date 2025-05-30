## README

for the model associated with:

Rulkov NF, Timofeev I, Bazhenov M (2004)
Oscillations in large-scale cortical networks: map-based model.
*J Comput Neurosci* **17**:203-23

**This model is written in C++.** The following instructions have worked under Redhat Linux 9 although the code should work with little or no modification on other platforms with ANSI compilers.

#### How to compile:

Change to the directory in which you have unzipped the archive file, then type:

```
g++ network2D.cpp -lm -O4 -ffast-math -o n.exe
```

#### How to run:

```
n.exe input2D.txt > tmp &
```

C++ code for simulating 2 layer 2D network (256x256 pyramidal cells - regular spiking type and 128x128 interneurons - fast spiking type) of map-based neurons.
With some small modifications this code was used to generate Fig 9A (right plot) in our JCNS paper.
Also there is a movie of spiral wave dynamics available (see pop-down window above model files in ModelDB).

#### Making the movie with matlab:

Matlab scripts which create images which are then used to create the movie:

```
import2D('graf_cx0')
movieMAP
```

To assemble the images into the movie and then save the movie:

```
for i=1:200; tmp=imread(['mov' num2str(i) '.png']); m(i)=im2frame(tmp); end;
movie2avi(m,'spiral2D.avi');
```

#### Additional note:

The simulations with 256x256 require > 1GB RAM (otherwise there is a huge time penalty for virtual memory use).


Correspondence:
Dr. Maxim Bazhenov, The Salk Institute,
10010 North Torrey Pines Road, La Jolla, CA 92037.
E-mail: [bazhenov@salk.edu](mailto:bazhenov@salk.edu?Subject=Novel%20model%20for%20study%20of%20large%20scale%20networks%202004)

---

2025-05-27 – Standardized to Markdown