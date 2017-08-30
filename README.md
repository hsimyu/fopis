# TDPIC: 三次元 EM Full Particle-in-cell code

## 概要
マルチグリッド法による宇宙機帯電解析を目指して開発中の三次元電磁粒子コードです。

## 必要ソフトウェアとライブラリ
- Cmake (ビルドツール)
    - windows では chocolatey でインストールするのが楽
- Boost
    - [boostjp](https://boostjp.github.io/) を参照
    - 環境変数 BOOST_ROOT にインストール場所を設定しておくこと
- Silo
    - 要HDF5とzlib
        - zlib は Cmake を利用することで Visual Studio でビルド可能
        - HDF5 も Cmake を利用することでビルド可能
            - [https://www.hdfgroup.org/downloads/hdf5/source-code/](https://www.hdfgroup.org/downloads/hdf5/source-code/)
        - Cmake用のビルドファイルを落としてきて、VS2017用のビルド設定を追加する
    - Siloを落としてきて, Copysilo.datを実行してsilo.h等を生成
        - 新しいMSVCでコンパイルする場合、snprintfのdefineがstdio.h内の定義と衝突するため、silo_win32_compatibility.hの15行目をコメントアウトする(4.10.2)
        - pdb_detect が perl を使って pdform.h を出力するため、perlがないとコンパイルが通らない
    - DevPropsみたいなのを設定する必要あり？
- MPIライブラリ
    - MSMPI (Windows), インストールすれば自動的に環境変数が設定されるはず
    - OpenMPI (Linux, macOS)
    - スパコン上では自動リンクされるため手動インストールは必要ない

- Doxygen (ドキュメント生成用, Optional)

## 初期ビルド

```
$ git clone https://snas1.rish.kuins.net:8787/gitlab/yamakawa-lab/tdpic.git ./tdpic
$ cd tdpic
$ git submodule init
$ git submodule update # googletestとpicojsonをgithubから取ってくる
$ mkdir build # cmakeでout-source ビルドするディレクトリ
$ cd build
$ cmake ..
```

もしくは
```
$ git clone https://snas1.rish.kuins.net:8787/gitlab/yamakawa-lab/tdpic.git ./tdpic
$ cd tdpic
$ ./scripts/setup.sh
```
でzlib, hdf5, siloのダウンロードとビルドまでをやります。

## 実行方法

## 物体定義について
ボクセルベースのobjファイルをサポート.
MagicaVoxelというフリーソフトで作るのが今のところ一番楽.
- MagicaVoxel @ ephtracy: https://ephtracy.github.io/

config/config.txtのobjファイル出力に関する部分を下記のように設定する.
- pivot = '0.5 0.5 0.0' -> '0.5 0.5 0.5'
- axis = 'XZ-Y' -> 'XYZ'
- optimize = '1' -> '0'

```
file_obj
(
    scale = '1 1 1'
    pivot = '0.5 0.5 0.5'
    tc_offset = '0.5 0.5'
    cw = '0'
    axis = 'XYZ'
    optimize = '0'
)
```
## ドキュメント
要 Doxygen

cmakeでmakefileを生成したあと、buildディレクトリで
```
$ make doc
```
すると生成されます。
```
$ open ./doc/html/index.html
```
すればブラウザでドキュメントが表示されます。
