# Running on SummitDev

### The makefile assumes that the kokkos source is in your home directory:

```
cd ~
git clone https://github.com/kokkos/kokkos
```

### Once that has been done load the correct environment
```
module load gcc
module load cuda
```
### And build
```
make
```

### Target Arch
To target a Kokkos execution space set `KOKKOS_DEVICES`
```
KOKKOS_DEVICES=Cuda
KOKKOS_DEVICES=OpenMP
KOKKOS_DEVICES=Serial
```
