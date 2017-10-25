CodeFEM - Codes for FEM

This depends on OpenCL eigen viennacl f2c

If you have no GPU, try amdapp-sdk(CPU's OpenCL).

* xmesh - A generator of FEM mesh.
  doc/xmesh.txt is how to use xmesh command.

* glirulus, glirulus() - Yet Another Matrix Applications' New Engine.
  doc/glirulus.txt is how to use glirulus command and glirulus() function.

* cl_bicgstab() - A linier solver using OpenCL.
  doc/cl_bicgstab.txt is how to use cl_bicgstab().

* psc98 - Parallel Software Contest '98 (PSC98)'s problems.
  doc/psc98.txt is how to use psc98.

* navierstokes() - A simulation of Navier-Stokes equations.
  doc/navierstokes.txt is how to use navierstokes().

### 1st. How to build all libraries.
$ cd bin/; ./mkall; cd ..
# mkdir -p /usr/include/est
# cp src//dataparallel.d/cl_bicgstab_kernel.cl /usr/include/est

### 2nd. How to build xmesh command.
$ cd src/xmesh.d; make; cd ../..

### 3rd. How to bulid glirulus command.
$ cd src/glirulus.d; make; cd ../..

### 4th. How to bulid psc98 command.
$ cd src/psc98.d; make; cd ../..
  
### 5th. How to build navierstokes sample.
# cd src/navierstokes.d; make; cd ../..



参考文献
A-Z
[ 1] JSPP PSC'98事務局: PSC98ホームページ, 並列処理シンポジウム JSPP, 1998.

[ 2] Richard Barrant/Michael Berry/Tony F. Chan/James Demmel/June Donato/Jack Dongarra/Victor Eijhout/Roldan Pozo/Charles Romine/Henk van der Vorst 著, 長谷川 里美/長谷川 秀彦/藤野 清次 訳: 反復法Templates, 朝倉書店, 1996.
あ
[ 3] 青木 尊之: OpenCL, 東京工業大学大学院情報理工学研究科GPUコンピューティングNo.14, 2013.

[ 4] 伊藤 祥司: 周期境界問題に対する並列計算機向きの前処理法(大規模線形方程式に対するスーパーコンピュータ向きの前処理法 第5章), 筑波大学工学研究科博士学位論文, 2001.

[ 5] 小国 力 編著: 行列計算ソフトウェア, 丸善株式会社, 1992.
か
[ 6] 金澤 靖: プログラムツールEigen 3次元コンピュータビジョン計算ハンドブック付録, 森北出版, 2016.

[ 7] 株式会社フィックスターズ: OpenCL入門, インプレスジャパン, 2010.

[ 8] 川原 睦人: 有限要素法流体解析, 日科技連, 1985.

[ 9] 河村 哲成: 二訂 流れのシミュレーションの基礎！, インデックス出版, 2007.

[10] 儀我 美一, 儀我 美保: 非線形偏微分方程式, 共立出版株式会社, 1999.

[11] 粂井康孝: 猫でもわかるC++プログラミング, ソフトバンク クリエイティブ株式会社, 2009.
さた
[12] 高見 利也: 流れ関数と渦度による二次元非圧縮流体の数値計算, 九州大学大学院システム情報科学府 仮想実験特論, 2010.

[13] 田口 圭一: 血液のキャビティ流れの有限要素解析, 高知工科大学工学部知能機械システム工学科卒業論文, 2004.

[14] 佃 良生: 複雑形状の領域における流れの数値計算について, 電気通信大学平成9年度情報工学科卒業論文, 1998.

[15] 佃 良生, 海津 聰: デローニー法によるFEM格子の自動形成, 数値流体力学会講演論文集, 1996.

[16] 薦田 登志矢, 三輪 忍, 中村 宏: OpenCLを用いたパイプライン並列プログラミングAPIの初期検討, 情報処理学会研究報告, 2011.

[17] 登坂 宣好, 角田 和彦: Navier-Stokes方程式の数値シミュレーション, 数理解析研究所講究録 第836巻, 1993.
な
[18] 中山 司: 流れ解析のための有限要素法入門, 東京大学出版会, 2008.

[19] 日本計算工学会 流れの有限要素法研究委員会 編: 続・有限要素法による流れのシミュレーション, シュプリンガー・フェアラーク東京, 2008.

[20] 日本計算工学会 流れの有限要素法研究委員会 編: 有限要素法による流れのシミュレーション, シュプリンガー・フェアラーク東京, 1998.
は
[21] 塙 与志夫, 小柳 義夫: 大規模疎行列係数連立一次方程式に対する前処理つき共役勾配法の並列化, ハイパフォーマンスコンピューティング 77-20,1999.
ま
[22] 村松 一弘, 鷲尾 巧, 土肥 俊: 並列マシンCenju2による非圧縮性流体解析, 情報処理学会研究報告ハイパフォーマンスコンピューティング Vol.1994, No.22, p.p9-16, 1994.
や
[23] 矢川 元基, 奥田 洋司, 中林 靖: 有限要素法流れ解析, 朝倉書店, 1998.
わ
[24] 渡部 善隆: 連立1次方程式の基礎知識, 数値解析チュートリアル2004資料, 2004.
--
ダウンロード
[ 1] については, https://web.kudpc.kyoto-u.ac.jp/doc/HPC-WG/PSC98/
[ 3] については, http://www.ocw.titech.ac.jp/index.php?module=General&action=DownLoad&file=20131226717065-480-0-52.pdf&type=cal&JWC=20131226717065
[ 4] については, https://tsukuba.repo.nii.ac.jp/?action=repository_action_common_download&item_id=8784&item_no=1&attribute_id=17&file_no=7
[ 6] については, https://www.morikita.co.jp/books/download/1574
[ 9] については, https://www.index-press.co.jp/books/science/ceslib2.pdf
[12] については, http://dogra.csis.oita-u.ac.jp/tkm/lecture/sim2010/pdf/streamFunction.pdf
[13] については, http://www.kochi-tech.ac.jp/library/ron/2003/2003mech/1040120.pdf
[14] については, http://kud.dip.jp/thesis.d/UEC_thesis.pdf
[14] の付属CD-ROMは, http://kud.dip.jp/thesis.d/UEC_thesis.CD_ROM/
[16] については, https://ipsj.ixsq.nii.ac.jp/ej/?action=repository_action_common_download&item_id=79275&item_no=1&attribute_id=1&file_no=1
[17] については, http://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/0836-13.pdf
[21] については, https://ipsj.ixsq.nii.ac.jp/ej/?action=repository_action_common_download&item_id=29555&item_no=1&attribute_id=1&file_no=1
[22] については, https://ipsj.ixsq.nii.ac.jp/ej/?action=repository_action_common_download&item_id=24237&item_no=1&attribute_id=1&file_no=1
[24] については, http://yebisu.cc.kyushu-u.ac.jp/~watanabe/RESERCH/MANUSCRIPT/TUTORIAL/leq.pdf

