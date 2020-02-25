data2gfs
*********
NCEP-GFSの初期値ファイル(sigma restart file)を作成するためのプログラム

概要
=====
ディレクトリ構成
-----------------
 ./exec 実行ファイル
 ./src  プログラムソース
 ./exp  実行シェルスクリプト

実効シェルスクリプト
---------------------
 ===========  ==============  ======= ========
 name         vertical level  format  user id
 ===========  ==============  ======= ========
 YOTC         hybrid level    grib1   userid=1
 GANAL        pressure level  grads   userid=2
 JRA25        pressure level  grib1   userid=3
 ERA-interim  hybrid level    grib1   userid=4
 ERA-interim  pressure level  grib1   userid=5
 ALERA2       sigma level     grads   userid=6
 ===========  ==============  ======= ========

対応フォーマット
-----------------
 grib1, grib2, grads形式

コンパイル
===========
NCEPのライブラリ(w3lib, bacio, g2lib, iplib, splib)が必要。

- http://www.nco.ncep.noaa.gov/pmb/docs/libs/
- http://www.nco.ncep.noaa.gov/pmb/codes/GRIB2/


更新履歴
=========
 2013.08.07 initial version 

Author Information
===================
tmiyachi 

- https://bitbucket.org/tmiyachi
- https://github.com/tmiyachi

