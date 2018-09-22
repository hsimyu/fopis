# FOPIS: Multigrid ES/EM Full Particle-in-cell Object-Plasma Interaction Simulator

## 概要
マルチグリッド法による宇宙機帯電解析を目指して開発中の三次元電磁粒子コードです。

## 必要ソフトウェアとライブラリ
- Cmake (ビルドツール)
    - windows では chocolatey でインストールするのが楽
- Boost
    - [boostjp](https://boostjp.github.io/) を参照
    - 環境変数 BOOST_ROOT にインストール場所を設定しておくこと
- HDF5とzlib
    - scripts/setup-on-local-unix.sh を利用するとよいです
    - WindowsにおけるParallel HDF5ビルドは少し大変
        - [https://www.hdfgroup.org/downloads/hdf5/source-code/](https://www.hdfgroup.org/downloads/hdf5/source-code/)
        - Cmake用のビルドファイルを落としてきて、VS2017用のビルド設定を追加する
- MPIライブラリ
    - Windowsの場合: MS MPIをインストールしてください
        - インストールすれば自動的に環境変数が設定されるはず
    - Unixの場合: OpenMPI (Linux, macOS)
        - brew install openmpi
    - スパコン上では自動リンクされるため手動インストールは必要ない

- Doxygen (ドキュメント生成用, Optional)

## Installation (UNIX)

### HDF5とzlibのインストールを自動でやる場合:

```bash
$ git clone https://snas1.rish.kuins.net:8787/gitlab/yamakawa-lab/tdpic.git ./tdpic
$ cd tdpic
```

インストール用のスクリプトがscripts/setup-on-local-unix.shにあります。

デフォルトでは$HOME/localにインストールを行います。
インストール先や、MPIコンパイラを変更したい場合はエディタで開いて先頭の5行の変数の値を変更してください。

```bash
$ ./scripts/setup-on-local-unix.sh
```

### HDF5とzlibのインストールが終わっている場合:

```bash
$ git clone https://snas1.rish.kuins.net:8787/gitlab/yamakawa-lab/tdpic.git ./tdpic
$ cd tdpic
$ git submodule init
$ git submodule update
$ mkdir build
$ cd build
$ cmake .. -DLOCAL_LIBRARYDIR=$HOME/local -DCMAKE_INSTALL_PREFIX=""
```

## ビルド & 実行方法

current directoryがbuildだと仮定します。

単一実行の場合:
```bash
$ make
$ ./tdpic
```
または
```bash
$ ./tdpic your_configuration_file.json
```
で実行可能です。
第一引数が与えられた場合、そのファイルを読み込み計算を行います。
第一引数が与えられない場合、同ディレクトリのinput.jsonを読み込み計算を行います。

別の場所へインストールする場合、make install時のDESTDIRで設定して下さい。

例:
```bash
$ make install DESTDIR=$HOME/Documents/simulation
```

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
