cd build
mkdir -p results/part1
git checkout 35b13008297e70b8ab60096bc4a7ac358b2ba9da
make -j `nproc`
# banana (~2k tris)
# building (~70k tris)

# beast (~60k tris)
# beetle (~60k tris)
# cow (~6k tris)
# maxplanck (~50k tris)
# peter (~40k tris)
# teapot (~2.5k tris)

# bench (~70k tris)
# blob (~200k tris)
# bunny (~30k tris)
# cbbunny (~30k tris)
# cbcoil (~7k tris)
# cbdragon (~100k tris)
# cbempty, cbgems are relatively trivial
# cblucy (~130k tris)
# cbspheres is relatively simple
# wall-e (~240k tris)

# plane / cube are trivial

# for part 1, these are the easiest ones
./pathtracer -r 480 360 -f results/part1/cube.png ../dae/simple/cube.dae
./pathtracer -r 480 360 -f results/part1/plane.png ../dae/simple/plane.dae
./pathtracer -r 480 360 -f results/part1/CBempty.png ../dae/sky/CBempty.dae
./pathtracer -r 480 360 -f results/part1/CBgems.png ../dae/sky/CBgems.dae
./pathtracer -r 480 360 -f results/part1/CBspheres.png ../dae/sky/CBspheres_lambertian.dae

mkdir -p results/part2

# for part 2, banana, cbcoil, cow, teapot are mildly challenging w/o taking way too long
# render bunny just to demonstrate that otherwise not feasible
if ! [ -f results/part2/timing ]; then
    printf "No BVH: \n" >> results/part2/timing 
    ./pathtracer -t 8 -r 480 360 -f results/part2/banana.png ../dae/keenan/banana.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cbcoil.png ../dae/sky/CBcoil.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cow.png ../dae/meshedit/cow.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/teapot.png ../dae/meshedit/teapot.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/bunny.png ../dae/sky/bunny.dae >> results/part2/timing
    
    git checkout f8f45db7cf1dff6892234f4bff6efcba71bcc06a
    make -j `nproc`
    printf "With BVH: \n" >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/banana.png ../dae/keenan/banana.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cbcoil.png ../dae/sky/CBcoil.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cow.png ../dae/meshedit/cow.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/teapot.png ../dae/meshedit/teapot.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/bunny.png ../dae/sky/bunny.dae >> results/part2/timing

    printf "With BVH: \n" >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/wall-e.png ../dae/sky/wall-e.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cbdragon.png ../dae/sky/CBdragon.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/cblucy.png ../dae/sky/CBlucy.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/blob.png ../dae/sky/blob.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/building.png ../dae/keenan/building.dae >> results/part2/timing
    ./pathtracer -t 8 -r 480 360 -f results/part2/maxplanck.png ../dae/meshedit/maxplanck.dae >> results/part2/timing
fi

mkdir -p results/part3

if ! [[ $* == *--skip3* ]]; then
    git checkout 89b6ccafdfeff28803df3b0fce8877db18c32129
    make -j `nproc`

    ./pathtracer -t 8 -s 64 -l 32 -m 6 -r 480 480 -f results/part3/dragon.png ../dae/sky/CBdragon.dae >> results/part3/dragon.timing
    ./pathtracer -t 8 -s 64 -l 32 -m 6 -H -r 480 480 -f results/part3/dragon_hemi.png ../dae/sky/CBdragon.dae >> results/part3/dragon_hemi.timing

    ./pathtracer -t 8 -s 64 -l 32 -m 6 -r 480 360 -f results/part3/coil.png ../dae/sky/CBcoil.dae >> results/part3/coil.timing
    ./pathtracer -t 8 -s 64 -l 32 -m 6 -H -r 480 360 -f results/part3/coil_hemi.png ../dae/sky/CBcoil.dae >> results/part3/coil_hemi.timing

    ./pathtracer -t 8 -s 64 -l 32 -m 6 -r 480 360 -f results/part3/bunny.png ../dae/sky/CBbunny.dae >> results/part3/bunny.timing
    ./pathtracer -t 8 -s 64 -l 32 -m 6 -H -r 480 360 -f results/part3/bunny_hemi.png ../dae/sky/CBbunny.dae >> results/part3/bunny_hemi.timing

    ./pathtracer -t 8 -s 1 -l 1 -m 6 -r 480 360 -f results/part3/bunny_1.png ../dae/sky/CBbunny.dae >> results/part3/bunny_1.timing
    ./pathtracer -t 8 -s 1 -l 4 -m 6 -r 480 360 -f results/part3/bunny_4.png ../dae/sky/CBbunny.dae >> results/part3/bunny_4.timing
    ./pathtracer -t 8 -s 1 -l 16 -m 6 -r 480 360 -f results/part3/bunny_16.png ../dae/sky/CBbunny.dae >> results/part3/bunny_16.timing
    ./pathtracer -t 8 -s 1 -l 64 -m 6 -r 480 360 -f results/part3/bunny_64.png ../dae/sky/CBbunny.dae >> results/part3/bunny_64.timing
fi

git checkout master
make -j `nproc`

mkdir -p results/part4
if ! [[ $* == *--skip4* ]]; then
    # direct vs indirect
    ./pathtracer -t 8 -s 128 -l 8 -m 1      -r 480 360 -f results/part4/empty_direct.png ../dae/sky/CBempty.dae >> results/part4/empty_direct.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6 -i 1 -r 480 360 -f results/part4/empty_indirect.png ../dae/sky/CBempty.dae >> results/part4/empty_indirect.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6      -r 480 360 -f results/part4/empty.png ../dae/sky/CBempty.dae >> results/part4/empty.timing
   
    ./pathtracer -t 8 -s 128 -l 8 -m 1      -r 480 360 -f results/part4/blob_direct.png ../dae/sky/blob.dae >> results/part4/blob_direct.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6 -i 1 -r 480 360 -f results/part4/blob_indirect.png ../dae/sky/blob.dae >> results/part4/blob_indirect.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6      -r 480 360 -f results/part4/blob.png ../dae/sky/blob.dae >> results/part4/blob.timing
   
    ./pathtracer -t 8 -s 128 -l 8 -m 1      -r 480 360 -f results/part4/spheres_direct.png ../dae/sky/CBspheres_lambertian.dae >> results/part4/spheres_direct.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6 -i 1 -r 480 360 -f results/part4/spheres_indirect.png ../dae/sky/CBspheres_lambertian.dae >> results/part4/spheres_indirect.timing
    ./pathtracer -t 8 -s 128 -l 8 -m 6      -r 480 360 -f results/part4/spheres.png ../dae/sky/CBspheres_lambertian.dae >> results/part4/spheres.timing
   
    # bunny bounce comparison
    ./pathtracer -t 8 -s 1024 -l 4 -m 0 -o 0 -r 480 360 -f results/part4/bunny_0_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_0_bounce.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 1 -o 0 -r 480 360 -f results/part4/bunny_1_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_1_bounce.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 2 -o 0 -r 480 360 -f results/part4/bunny_2_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_2_bounce.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 3 -o 0 -r 480 360 -f results/part4/bunny_3_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_3_bounce.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 4 -o 0 -r 480 360 -f results/part4/bunny_4_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_4_bounce.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 5 -o 0 -r 480 360 -f results/part4/bunny_5_bounce.png ../dae/sky/CBbunny.dae >> results/part4/bunny_5_bounce.timing

    ./pathtracer -t 8 -s 1024 -l 4 -m 0 -r 480 360 -f results/part4/bunny_0_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_0_accum.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 1 -r 480 360 -f results/part4/bunny_1_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_1_accum.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 2 -r 480 360 -f results/part4/bunny_2_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_2_accum.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 3 -r 480 360 -f results/part4/bunny_3_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_3_accum.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 4 -r 480 360 -f results/part4/bunny_4_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_4_accum.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 5 -r 480 360 -f results/part4/bunny_5_accum.png ../dae/sky/CBbunny.dae >> results/part4/bunny_5_accum.timing

    # roulette
    ./pathtracer -t 8 -s 1024 -l 4 -m 0   -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_0_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_0_roulette.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 1   -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_1_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_1_roulette.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 2   -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_2_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_2_roulette.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 3   -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_3_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_3_roulette.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 4   -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_4_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_4_roulette.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 100 -z 1 -x 0.65 -r 480 360 -f results/part4/bunny_5_roulette.png ../dae/sky/CBbunny.dae >> results/part4/bunny_5_roulette.timing

    # sample rate
    ./pathtracer -t 8 -s 1    -l 4 -m 5 -r 480 360 -f results/part4/bunny_1_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_1_sample.timing
    ./pathtracer -t 8 -s 2    -l 4 -m 5 -r 480 360 -f results/part4/bunny_2_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_2_sample.timing
    ./pathtracer -t 8 -s 4    -l 4 -m 5 -r 480 360 -f results/part4/bunny_4_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_4_sample.timing
    ./pathtracer -t 8 -s 8    -l 4 -m 5 -r 480 360 -f results/part4/bunny_8_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_8_sample.timing
    ./pathtracer -t 8 -s 16   -l 4 -m 5 -r 480 360 -f results/part4/bunny_16_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_16_sample.timing
    ./pathtracer -t 8 -s 64   -l 4 -m 5 -r 480 360 -f results/part4/bunny_64_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_64_sample.timing
    ./pathtracer -t 8 -s 1024 -l 4 -m 5 -r 480 360 -f results/part4/bunny_1024_sample.png ../dae/sky/CBbunny.dae >> results/part4/bunny_1024_sample.timing
fi

mkdir -p results/part5
if ! [[ $* == *--skip5* ]]; then
    ./pathtracer -t 8 -s 2048 -l 1 -m 5 -v 1 -r 480 360 -f results/part5/bunny_adaptive.png ../dae/sky/CBbunny.dae >> results/part5/bunny_adaptive.timing
    ./pathtracer -t 8 -s 2048 -l 1 -m 5 -v 1 -r 480 360 -f results/part5/spheres_adaptive.png ../dae/sky/CBspheres_lambertian.dae >> results/part5/spheres_adaptive.timing
fi