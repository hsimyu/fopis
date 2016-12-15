# TDPIC: 三次元 EM Full Particle-in-cell code

## 概要
マルチグリッド法による宇宙機帯電解析を目指して開発中の三次元電磁粒子コードです。

## 使い方
ビルドツール: Cmake
ドキュメント生成ツール: Doxygen

`
$ git clone https://snas1.rish.kuins.net:8787/gitlab/yamakawa-lab/tdpic.git ./tdpic
$ cd tdpic
$ git submodule init
$ git submodule update # googletestとpicojsonをgithubから取ってくる
$ mkdir build # cmakeでout-source ビルドするディレクトリ
$ cd build
$ cmake ..
`

## ドキュメント
要 Doxygen

cmakeでmakefileを生成したあと、buildディレクトリで
`$ make doc`
すると生成されます。
`$ open ./doc/html/index.html`
すればブラウザでドキュメントが表示されます。
