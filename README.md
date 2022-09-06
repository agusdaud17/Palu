# Palu
Analisis kapasitas dukung vertikal pondasi tiang (bor dan pancang)

## Fitur

1. Kapasitas dukung vertikal pondasi
2. Estimasi penurunan tiang dengan kurva $t-z$ (Chen and Kulhawy, 2002).
3. Analisis kapasitas multi-method

## Metode: Pondasi Bor (Drilled Pile)
### Pasir (Cohesionless soil)

1. Tahanan selimut, $f_s$ (*skin friction*)
    - Reese and O'Neil
    - Meyerhof
    - Reese and Wright
    - Touma and Reese
    - Quiros and Reese
    - FHWA
    - Eurocode
    - API
    
2. Tahanan ujung, $f_b$ (*tip resistance*)
    - Meyerhof
    - Reese and O'Neil
    - Reese and Wright
    - Touma and Reese (N,dp)
    - API
    - Eurocode

### Lempung/Lanau Plastis (Cohesive soil)

1. Tahanan selimut, $f_s$ (*skin friction*)
    - Skempton
    - Reese and O'Neil
    - FHWA
    - API
    
2. Tahanan ujung, $f_b$ (*tip resistance*)
    - Skempton
    - Reese and O'Neil

**Service limit state:** using normalized load-deformation curves $(t-z)$, Chen and Kulhawy (2002)

## Metode: Pondasi Tiang Pancang (*Driven Pile*)

### Pasir (Cohesionless soil)

1. Tahanan selimut, $f_s$ (*skin friction*): Meyerhof
    
2. Tahanan ujung, $f_b$ (*tip resistance*) : Meyerhof

### Lempung/Lanau Plastis (Cohesive soil)

1. Tahanan selimut, $f_s$ (*skin friction*): Skempton
    
2. Tahanan ujung, $f_b$ (*tip resistance*) : API

## Input data
Input data yang diperlukan, berupa format excel. Dapat dilihat pada file studi_kasus.xlsx

## IDE/Text Editor
Untuk menggunakan module ini direkomendasikan menggunakan Jupyter Notebook/ Jupyter Lab.

# Pondasi tiang bor (Drilled Pile)

## Single method

```
from Palu import Drilled

p1 = Drilled(
    Lp=43, # kedalaman dasar tiang
    D=300, # diameter pondasi (cm)
    GWL=15.3, # muka air tanah (m)
    ignore_Ltop=1.5, # panjang (m) tahanan friksi yang diabaikan pada bagian atas tiang
    ignore_Lbot=3,# panjang (m) tahanan friksi yang diabaikan pada bagian bawah tiang
    dz=0.5, # panjang element pias/tinjauan
)

p1.read_excel('studi_kasus.xlsx') # input data

p1.solve(
    FS=2.5, # faktor aman
    Qb = ["Reese and Wright","Skempton"], # metode analisis tahanan ujung [tanah pasir, tanah kohesif]
    Qs = ["Meyerhof","Reese and O'Neil"]# metode analisis tahanan selimut [tanah pasir, tanah kohesif]
)
```

Menampilkan hasil:
```
p1.info()
```
```
Qs = 22405.88 kN
Qb = 21300.98 kN
Wp = 5371 kN
Qult = 43706.86 kN
Qa = 17482.74 kN
Q uplift = 8333.0 kN
```

### Kurva $t-z$

```
# Data hasil test pembebanan pondasi
data_P = [0,17800,26660,35590,40033,44500]
data_S = [0,0.254,0.762,5.8,10.7,19.0]
data = [data_P,data_S]

# Kurva t-z
p1.TZ_curve(S=20,test_data=data)
```

![Kurva t-z](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgXmRRyos4glH_bJk4y48wX9dZentCxhMECH2FusOFhjs1VjKnm4jgN1e2CxJK7b73m4S9RRs9SgXBoW8JYUI0qbQxaPNG9UexLR9Kq6FL5IBOIMHxwkstgk_I6mMlTcSK1eVBCLtsplA_U5yjKxOPJo1drniaEQJAv9IL6GNxELBjUehMdD8kY_pa-UQ/s600/kurva_t_z.png)

## Multi-method

```
from Palu import Drilled

pm = Drilled(Lp=43,D=300,GWL=15.3,ignore_Ltop=1.5,ignore_Lbot=3,dz=0.5)
pm.read_excel('studi_kasus.xlsx')
pm.Msolve(FS=2.5,ignore_fs_end=False)
pm.plot_result(param="Qa")
```
Hasil:

![Result](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhz6rMQWsDITYE7UTpjb_fRuHIL4TjcWRazWuzg0bKcWviGClqFktU_VtTAKqbhbKXLY6D5PqypuTM6sTiqIOxqye_LWEceq6_spK3ySb8MrA1UfOSb_FT1wMNvok7hdr7SZJp8386B6AM66-4SzZv5BMjpQFGgf5OT2_xXU5gaF28UrEbWd9eTDxQKPA/w640-h534/multi_method_drilled.png)

# Pondasi tiang pancang (Driven Pile)
Contoh:
```
from Palu import Driven

pd = Driven(Lp=43,D=300,GWL=15.3,ignore_Ltop=1.5,ignore_Lbot=3,dz=0.5)
pd.read_excel('studi_kasus.xlsx')
pd.solve(FS=2.5,ignore_fs_end=False)

pd.info()
```
Hasil:
```
Qs = 34949.95 kN
Qb = 21300.98 kN
Wp = 5371 kN
Qult = 56250.93 kN
Qa = 22500.37 kN
Q uplift = 12096.0 kN
```
