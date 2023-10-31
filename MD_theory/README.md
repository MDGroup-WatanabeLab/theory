# Molecular Dynamics 理論
&emsp;本ドキュメントでは、LAMMPSやVASPを用いるにあたり、必要となるであろうMDの理論を説明する。参考文献は[先生が書いた本](https://www.morikita.co.jp/books/mid/092251)と筆者が作った[真空ゼミのスライド](https://waseda.app.box.com/file/1080558360236)である。後者については間違っている部分もあるので、前者を主に参考に

## 目次
1. [概要](#1-概要)
2. [運動方程式](#2-運動方程式)
3. [ポテンシャル](#3-ポテンシャル)
4. [系と周期的境界条件](#4-系と周期的境界条件)
5. [構造最適化](#5-構造最適化)
6. [温度制御](#6-温度制御)
7. [圧力制御](#7-圧力制御)

※クリックで各セクションへ移動可能

## 1. 概要
&emsp;MDとは、Molecular Dynamicsの頭文字をとっており、__分子動力学法__ と訳される。古典力学を用いたシミュレーション方法であり、主に原子間に働く力をポテンシャルから算出し、原子を動かしながら物性を予測するシミュレーションである。昨今では、今まで使われていた経験的ポテンシャルに代わり、機械学習ポテンシャルを用いた計算が主流となってきている。本ドキュメントでは、機械学習ポテンシャルには触れず、MDシミュレーションがどのような流れで行われているかについて述べる。  
　まずは、力学内でのMDの位置づけを見ていく。下の図を見てほしい。  
<div align = "center">
<img width="681" alt="スクリーンショット 2023-10-25 095231" src="https://github.com/MDGroup-WatanabeLab/image_for_mdpython/assets/138444525/b6925c81-6785-4595-aca7-ae0396bcd96a"> 
</div>   
この図は有名だが、量子力学や古典力学、連続体力学のシミュレーションの種類と、時間スケール、扱える原子数を示している。MDであれば、1ps~1nsスケールのシミュレーションを扱えるという意味である。 縦軸は計算コストではないので注意が必要である。また、研究の観点からの位置づけも見ていく。  
<div align = "center">
<img width="402" alt="スクリーンショット 2023-10-25 101730" src="https://github.com/MDGroup-WatanabeLab/image_for_mdpython/assets/138444525/08a48783-18cc-435b-acdd-664f77243190">
</div>  
このように、実験を行う前に理論に基づいた計算を行うことで、実験コストを最小限に抑えたり、危険性や毒性がある物質の取り扱いが容易となるのである。さて、本編に入る前に、MDシミュレーションの大まかな流れを最後に述べる。  

```mermaid
graph TD;
    原子の配置とポテンシャルの設定 --> 原子間に働く力を計算;
    原子間に働く力を計算 --> 運動方程式を解く;
    運動方程式を解く --> 原子を動かす;
    原子を動かす --> 系のエネルギーの計算;
    原子を動かす --> 原子間に働く力を計算;
```  

このフローチャートの中身についてこれから説明していくこととする。

## 2. 運動方程式
&emsp;順番は前後するが、まずは運動方程式について説明する。運動方程式がなぜ必要なのかというと、MDでは原子の座標の情報が必要不可欠だからである。原子にかかる力がわかれば、加速度がわかり、速度や $\Delta t$ 秒後の位置も求まるというのは想像に容易いだろう。  
<div align = "center">
<img width="609" alt="スクリーンショット 2023-10-25 113039" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/73f5f0fb-452a-4045-8459-089fe2c09813">
</div>  
　さて、ポテンシャルから力を計算をすると述べた。大学物理の復習となるが、式で表すと次のようになる。  

$$  
\boldsymbol{F} = m \frac{d^2 \boldsymbol{r}}{dt^2}=-\frac{\partial V}{\partial \boldsymbol{r}}
$$  

ポテンシャルエネルギーは $V$ とした。このように、運動方程式を解くにはポテンシャルエネルギーが必須なのである。では、この運動方程式はどのようにMDシミュレーションで使われ、解かれるのか？  

### 2.1 Verlet法
　代表的な方法として、Verlet法と速度Verlet法がある。まずは、Verlet法において用いられる $\Delta t$ 秒後の位置を表す式を示す。  

$$  
\boldsymbol{r_i}(t+\Delta t) = 2\boldsymbol{r_i}(t)-\boldsymbol{r_i}(t-\Delta t)+(\Delta t)^2 \frac{\boldsymbol{F_i}(t)}{m} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.1.1)
$$  

太字はベクトル、添え字の $i$ は $i$ 番目の原子という意味である。つまり、(2.1.1)式を次のように繰り返し用いれば、原子座標を求めることが可能である。  

```mermaid
graph TD;
    初期値としてt秒後とt-Δt秒後の位置を与える --> 2.1.1式を用いてt+Δt秒後の位置を求める;
    2.1.1式を用いてt+Δt秒後の位置を求める --> 計算終了;
    2.1.1式を用いてt+Δt秒後の位置を求める --> 初期値としてt秒後とt-Δt秒後の位置を与える;
```  

では、この式はどのように導かれるのか？  
　t+Δt秒後の位置とt-Δt秒後の位置をそれぞれΔt=0の周りで2次の項までテイラー展開すると、

$$  
\boldsymbol{r_i}(t+\Delta t)=\boldsymbol{r_i}+\Delta t\cdot
\frac{d\boldsymbol{r_i}}{dt} +(+\Delta t)^2 \frac{1}{2}\frac{d^2\boldsymbol{r_i}}{dt^2} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.1.2)
$$  

$$  
\boldsymbol{r_i}(t-\Delta t)=\boldsymbol{r_i}-\Delta t\cdot
\frac{d\boldsymbol{r_i}}{dt} +(-\Delta t)^2 \frac{1}{2}\frac{d^2\boldsymbol{r_i}}{dt^2}+ O(\Delta t^3)\ \ \ \ \ \ \ \ (2.1.3)
$$  

(2.1.2)式と(2.1.3)式を加算すると、  

$$  
\boldsymbol{r_i}(t+\Delta t)+\boldsymbol{r_i}(t-\Delta t)=2\boldsymbol{r_i}+(+\Delta t)^2\cdot\frac{d^2\boldsymbol{r_i}}{dt^2}+O(\Delta t^3)\ \ \ \ \ (2.1.4)
$$  

変形すると、  

$$  
\frac{d^2\boldsymbol{r_i}}{dt^2}=\frac{\boldsymbol{r_i}(t+\Delta t)+\boldsymbol{r_i}(t-\Delta t)-2\boldsymbol{r_i}}{(\Delta t)^2}++O(\Delta t)\ \ \ \ \ \ \ \ \ (2.1.5)
$$  

(2.1.5)式をNewtonの運動方程式に代入すると、  

$$  
\frac{\boldsymbol{F_i}(t)}{m}=\frac{\boldsymbol{r_i}(t+\Delta t)+\boldsymbol{r_i}(t-\Delta t)-2\boldsymbol{r_i}}{(\Delta t)^2}++O(\Delta t)\ \ \ \ \ \ \ \ \ (2.1.6)
$$  

(2.1.6)式を変形すれば、  

$$  
\boldsymbol{r_i}(t+\Delta t) = 2\boldsymbol{r_i}(t)-\boldsymbol{r_i}(t-\Delta t)+(\Delta t)^2 \frac{\boldsymbol{F_i}(t)}{m} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.1.1)
$$  

と導かれる。

### 2.2 速度Verlet法
&emsp;速度Verlet法はVerlet法の発展形と思ってほしい。式の形は、  

$$  
\boldsymbol{r_i}(t+\Delta t) = 2\boldsymbol{r_i}(t)-\boldsymbol{r_i}(t-\Delta t)+(\Delta t)^2 \frac{\boldsymbol{F_i}(t)}{m} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.1.1)
$$  

$$  
\boldsymbol{v_i}(t+\Delta t) = \boldsymbol{v_i}(t)+\Delta t\cdot \frac{\boldsymbol{F_i}(t)}{m}+\frac{\Delta t}{2m}\left[ \boldsymbol{F_i}(t+\Delta t)-\boldsymbol{F_i}(t) \right]+O(\Delta t^3)\ \ \ \ \ \ \ \ (2.2.1)
$$  

新たに、(2.2.1)式が追加されている。Verlet法と異なり、速度の情報も導かれることがわかる。この2式を使う流れとしては、  

```mermaid
graph TD;
    初期値としてt秒後の位置と速度を与える --> 2.1.1式を用いてt+Δt秒後の位置を求める;
    2.1.1式を用いてt+Δt秒後の位置を求める --> 2.2.1式を用いてt+Δt秒後の速度を求める;
    2.2.1式を用いてt+Δt秒後の速度を求める --> 初期値としてt秒後の位置と速度を与える;
```  

　
では、式の導出を行う。その際、__前進差分近似と後進差分近似__ という数学の手法を知る必要があるので、先に解説する。t+Δt秒後の位置とt-Δt秒後の位置をそれぞれΔt=0の周りで1次の項までテイラー展開すると、

$$  
\boldsymbol{r_i}(t+\Delta t)=\boldsymbol{r_i}+\Delta t\cdot
\frac{d\boldsymbol{r_i}(t)}{dt} + O(\Delta t^2)\ \ \ \ \ \ \ \ (2.2.2)
$$  

$$  
\boldsymbol{r_i}(t-\Delta t)=\boldsymbol{r_i}-\Delta t\cdot
\frac{d\boldsymbol{r_i}(t)}{dt} + O(\Delta t^2)\ \ \ \ \ \ \ \ (2.2.3)
$$  

となる。それぞれの式を変形すると、  

$$  
\frac{d\boldsymbol{r_i}}{dt} = \frac{\boldsymbol{r_i}(t+\Delta t)-\boldsymbol{r_i}(t)}{\Delta t} + O(\Delta t)\ \ \ \ \ \ \ \ (2.2.4)
$$  

$$  
\frac{d\boldsymbol{r_i}}{dt} = \frac{\boldsymbol{r_i}(t)-\boldsymbol{r_i}(t-\Delta t)}{\Delta t} + O(\Delta t)\ \ \ \ \ \ \ \ (2.2.5)
$$  

(2.2.4)式は後ほど使用する。では、まず、t+Δt秒後の位置と速度をΔt=0の周りで2次までテイラー展開すると、  

$$  
\boldsymbol{r_i}(t+\Delta t)=\boldsymbol{r_i}+\Delta t\cdot\boldsymbol{v_i}(t) +(+\Delta t)^2 \frac{\boldsymbol{F_i}(t)}{m} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.2.6)
$$  

$$  
\boldsymbol{v_i}(t+\Delta t)=\boldsymbol{r_i}+\Delta t\cdot
\frac{\boldsymbol{F_i}(t)}{m} +(+\Delta t)^2 \frac{1}{2m}\frac{d\boldsymbol{F_i}}{dt} + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.2.7)
$$  

(2.2.4)式を(2.2.7)式の右辺第3項に適用すると、  

$$  
\boldsymbol{v_i}(t+\Delta t)=\boldsymbol{r_i}+\Delta t\cdot
\frac{\boldsymbol{F_i}(t)}{m} +\frac{\Delta t}{2m}\left[\boldsymbol{F_i}(t+\Delta t)-\boldsymbol{F_i}(t)\right] + O(\Delta t^3)\ \ \ \ \ \ \ \ (2.2.1)
$$  

これで導出は完了である。

## 3. ポテンシャル
&emsp;前の節で述べた通り、運動方程式の利用にはポテンシャルエネルギーが必要不可欠である。昨今では、機械学習ポテンシャルを作成し、計算を行う場合も出てきているが、それまでは経験的ポテンシャルというものが用いられていた。検索するだけで様々な種類のポテンシャルが出てくるが、本ドキュメントでは、いくつか抜粋して紹介するにとどめる。なお、機械学習ポテンシャルについては別のドキュメントで触れることとする。  
　紹介の前に、ポテンシャルそのものについて説明する。頭の中で、原子が64個ある系を想像してほしい。__原子1個に注目したとき、__ 相互作用はどうなっているだろう？2体系は63個あり、3体系は1935個ある。4体系等も考えていくと大変なことになる。  

```math
V_{total} = \sum_{j=2}^N \sum^{j-1}_{i=1}V_2(\boldsymbol{r_i, r_j})+\sum_{k=3}^N\sum_{j=2}^{k-1} \sum^{j-1}_{i=1}V_3(\boldsymbol{r_i, r_j, r_k}) + \cdots + V_N(\boldsymbol{r_1, r_2,\cdots, r_N})
``` 

相互作用は系の拡張によりどんどん大きくなっていくので、近似を行うことで系の拡張にも対応できるようにしなければならない。詳しい近似の内容は、ポテンシャルや適用したい系などにより異なるので適宜確認しなければならない。

### 3.1 Coulomb ポテンシャル
&emsp;電荷をもつ物体同士に働く相互作用であり、  

$$  
V=\frac{1}{4\pi\epsilon}\frac{Z_1Z_2}{|\boldsymbol{r_{21}}|}\ \ \ \ \ \ \ \ (3.1.1)
$$  

$\epsilon$ は誘電率、$\boldsymbol{r_{21}}$ は粒子間距離、Zは電荷を表す。電磁気学で学習済みだと思うので特に詳しい説明はしない。

### 3.2 Stilinger-Weber ポテンシャル
&emsp;SWポテンシャルと略されることもあるポテンシャルで、主にダイアモンド構造を取るⅣ族原子やその混晶、共有結合結晶に適用される。式としては、  

$$  
V(r)=\epsilon\sum_i\sum_{j>i}f_2\left(\frac{r_{ij}}{\sigma}\right)+\epsilon\sum_i\sum_{j>i}\sum_{k>j}f_3\left(\frac{\boldsymbol{r_i}}{\sigma}, \frac{\boldsymbol{r_j}}{\sigma}, \frac{\boldsymbol{r_k}}{\sigma}\right)\ \ \ \ \ \ \ \ (3.2.1)
$$  

$$  
f_2(r)=A\left(Br^{-p}-r^{-q}\right) exp \left( \frac{\delta}{r-a} \right)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3.2.2)
$$  

$$  
f_3\left(\boldsymbol{r_i}, \boldsymbol{r_j}, \boldsymbol{r_k}\right)= h(r_{ij}, r_{ik}, \theta_{jik})+h(r_{ji}, r_{jk}, \theta_{ijk})+h(r_{ki}, r_{kj}, \theta_{ikj})\ \ \ \ (3.2.3)
$$  

$$  
h(r_{ij}, r_{ik}, \theta_{jik})=\lambda exp\left[ \left(\frac{\gamma}{r_{ij}-a}\right)+\left(\frac{\gamma}{r_{ik}-a}\right) \right]\left( cos\theta_{jik}+\frac{1}{3} \right)^2\ \ \ \ (3.2.4)
$$  

各パラメータについては以下の表のとおり。 
|記号|単位|意味|
|:--:|:--:|:--:|
|$a$|Å|カットオフ距離|
|$\epsilon$|eV|結合の強さによるパラメータ|
|$\sigma$|Å|原子の大きさに依るパラメータ|
|$A, B, p, q, \delta, \lambda, \gamma$|$L^0M^0T^0$|フィッティングで決定されるパラメータ|

2体系と3体系の２つで構成されているため、4体系以降は考えない。  

### 3.3 Lennard-Jones ポテンシャル
&emsp;アルゴンなどの希ガスを対象とした2体間相互作用ポテンシャルであり、式は、

$$  
V(r)=4\epsilon\left[ \left( \frac{\sigma}{r} \right)^6 - \left( \frac{\sigma}{r} \right)^{12}\right] \ \ \ \ \ \ \ \ (3.3.1)
$$  

各パラメータについては以下の表のとおり。 
|記号|単位|意味|
|:--:|:--:|:--:|
|$\epsilon$|eV|エネルギーの最小値の値|
|$\sigma$|Å|原子の大きさに依るパラメータ|

ポテンシャルのカーブは次のようになる。  
<div align = "center">
<img width="286" alt="スクリーンショット 2023-10-25 141014" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/43f1cea6-f105-4630-a9d6-c88e29483321">
</div>

### 3.4 Born-Mayer-Huggins ポテンシャル
&emsp;イオン結晶を対象としたポテンシャルであり、式は、

$$  
V(r)= \frac{Z_iZ_j e^2}{4\pi\epsilon r}+A_{ij}b\cdot exp\left( \frac{\sigma_i+\sigma_j-r_{ij}}{\rho} \right)-\frac{C_{ij}}{r_{ij}^6}-\frac{D_{ij}}{r_{ij}^8} \ \ \ \ \ \ (3.4.1)
$$  

各パラメータについては以下の表のとおり。 
|記号|単位|意味|
|:--:|:--:|:--:|
|$e$|$C$|電気素量(1.60217662)|
|$Z_i, Z_j$||イオンの価数|
|$r_{ij}$|Å|イオン間距離|
|$A_{ij}$||ポーリング因子|
|$b$||反発力に依るパラメータ|
|$\sigma_i, \sigma_j$||原子の大きさに依るパラメータ|
|$\rho$||ソフトネスパラメータ|

## 4. 系と周期的境界条件
&emsp;では、計算の対象となる系にはどのような制限や条件が課されているか見てみよう。まず、実際に使用する構造の図を示す。  
<div align = "center">  
<img width="543" alt="スクリーンショット 2023-10-19 170018" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/72b9a411-267b-4286-ba3d-4fa6f4ccb98a">
</div>  
これだけ見ると、「え？こんだけの原子数の構造を計算して何になるの？意味なくね？」と思うかもしれない。しかし、計算を行う際には、この構造に対し、周期的境界条件というものが課される。  
　この節では、「既約表現論」の内容を復習しつつ、MD計算と対象となる系、その後に周期的境界条件について説明する。  

### 4.1 系
&emsp;さて、既約表現論や材料科学概論の内容を復習していこう。MD計算では、単体や多元系のアモルファスや結晶のエネルギーや力などの物理量を計算する。計算班では、Materials Projectや[自作のプログラム群](https://github.com/MDGroup-WatanabeLab/mdpython)で結晶構造を作成している。なので、三方晶や立方晶、六方晶などの知識が必要となってくるのである。したがって、この節では、7晶系やfccなどのよく出てくる結晶系の構成について簡単に復習する。  
　ではまず、7晶系である。下の表を見ながら思い出してほしい。


|晶系|軸|軸角|ブラベー格子|
|:-:|:-:|:-:|:-:|
|三斜|$a\not ={b\not ={c}}$|$\alpha\not ={\beta\not ={\gamma}}$|P|
|単斜|$a\not ={b\not ={c}}$|$\alpha=\gamma=90\degree$|P, C|　
|直方（斜方）|$a\not ={b\not ={c}}$|$\alpha=\beta=\gamma=90\degree$|P, C, I, F|　
|正方|$a=b$|$\alpha=\beta=\gamma=90\degree$|P, I|　
|三方|$a=b=c$|$\alpha=\beta=\gamma$|R|　
|六方|$a=b$|$\alpha=\beta=90\degree, \gamma=120\degree$|P|　
|立方（等軸）|$a=b=c$|$\alpha=\beta=\gamma=90\degree$|P, I, F|　

※ P：単純、C：底心、I：体心、F：面心、R：菱面体 を表す  

次に、代表的な結晶構造について見ていく。

|構造名|原子配置|例|
|:-:|:-:|:-:|
|NaCl型|$A(0, 0, 0), A(0, \frac{1}{2}, \frac{1}{2}), A(\frac{1}{2}, 0, \frac{1}{2}), A(\frac{1}{2}, \frac{1}{2}, 0)$<br> $B(\frac{1}{2}, 0, 0), B(0, \frac{1}{2}, 0), B(0, 0, \frac{1}{2}), B(\frac{1}{2}, \frac{1}{2}, \frac{1}{2})$|NaCl, KCl|
|FCC（面心立方）|$(0, 0, 0), (\frac{1}{2}, \frac{1}{2}, 0), (\frac{1}{2}, 0, \frac{1}{2}), (0, \frac{1}{2}, \frac{1}{2})$|Al, Cu|
|BCC（体心立方）|$(0, 0, 0), (\frac{1}{2}, \frac{1}{2}, \frac{1}{2})$|Fe, Cr|
|HCP（六方最密充填構造）|$(0, 0, 0), (\frac{1}{3}, \frac{2}{3}, \frac{1}{2})$|Mg, Zr|
|ダイアモンド型|$(0, 0, 0), (0, \frac{1}{2}, \frac{1}{2}), (\frac{1}{2}, 0, \frac{1}{2}), (\frac{1}{2}, \frac{1}{2}, 0)$<br>$(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}), (\frac{1}{4}, \frac{3}{4}, \frac{3}{4}), (\frac{3}{4}, \frac{1}{4}, \frac{3}{4}), (\frac{3}{4}, \frac{3}{4}, \frac{1}{4})$|C, Si|
|CsCl型|$A(0, 0, 0), B(\frac{1}{2}, \frac{1}{2}, \frac{1}{2})$|CsCl, CsBr|

このような構造はわかっておいた方が良いだろう。計算する際は、上記の単位胞をxyz方向にそれぞれ2倍、3倍して拡張した系を使用する。実際、先ほど示した結晶はGeのダイアモンド型構造を2×2×2したものである。

### 4.2 周期的境界条件
&emsp;では、本題の周期的境界条件について説明する。この周期的境界条件を用いる一番の理由は、__表面効果を無視できる__ からである。表面効果というのは、文字通り、結晶などの表面においてのみ見られる特別な現象のことをいう。SiやGeの結晶構造であれば、表面付近の原子間の結合数は4個から1個や2個に減ることとなり、通常よりも原子が動きやすくなってしまう。それは、バルクの計算としては不適切な結果が得られてしまうと考えられる。  
　したがって、MD計算を行う際は、原子が格子を飛び出した場合はその逆から原子が入ってくるといった形を取ることで、格子が無限に続いている状態を再現しているのである。下の図の赤い枠の格子をオリジナルセル、その他の青枠の格子をイメージセルと呼ぶ。  

<div align = "center">  
<img width="654" alt="スクリーンショット 2023-10-30 161649" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/d8a070d8-9bd8-4fd6-a3ed-efda9fc355cf">
</div>  

&emsp;しかし、この話を聞いて勘のいい人は、「え？無限に格子が続いていると仮定したら、エネルギーの値やばくね？」と考えるだろう。そのとおりであり、このままではエネルギーの値が非常に大きくなり、計算コストも爆上がりしてしまうことはすぐにわかる。この問題を解決するために用いられるのが __カットオフ距離とEwald法__ である。  

### 4.3 エネルギーのカットオフ
&emsp;というわけで、まずはカットオフについて説明する。カットオフ距離については、少々注意する必要がある。下の図を見てほしい。

<div align = "center">  
<img width="628" alt="スクリーンショット 2023-10-30 171948" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/62a63227-0942-41c3-be1c-1e6d0006efd7">
</div>  

このように、カットオフ距離が大きすぎると、ポテンシャルエネルギーが重複されて計算されてしまうため、単位胞の大きさの半分以下の長さにすることが好ましいと言える。ここで、勘のいいひとは、「え？そんなことしたら逆にエネルギーが本来よりも小さく出てくるのでは？」と考えるだろう。これを解決するのが、__Ewald法__ である。  

### 4.4 Ewald法
&emsp;というわけで、エネルギーのカットオフにおいて必要不可欠なEwald法について説明する。先述したクーロンポテンシャルの式は、

$$  
V=\frac{1}{4\pi\epsilon}\frac{Z_1Z_2}{|\boldsymbol{r_{21}}|}\ \ \ \ \ \ \ \ (3.1.1)
$$  

であった。この式から、クーロンポテンシャルは距離の逆数 $\frac{1}{r}$ に比例した関数であるとわかる。しかし、この関数は下記の通りの形をとり、減衰しづらいことがわかる。  

<div align = "center">  
<img width="463" alt="スクリーンショット 2023-10-30 181318" src="https://github.com/MDGroup-WatanabeLab/theory/assets/138444525/5efe4134-d392-461d-bf39-3a8aea0c2ab1">
</div>  

減衰しづらいということは、カットオフ距離より外側でエネルギーを０と近似したときの誤差が大きくなってしまうと言える。したがって、クーロンポテンシャルのような減衰しづらいポテンシャルは計算方法を工夫しなければならないとわかる。この工夫こそがEwald法である。今後、この節では、全ポテンシャルエネルギーを $V(r)$ として説明していく。  
　今回は、電荷中性条件が成り立っており、N個のイオンが存在する体積 $L^3$ のユニットセルを考える。まず、前述の通り、ポテンシャルは減衰しやすいものと減衰しづらいものの2つに分けられるので、次の式のように分解する。  

$$  
V(r) = V^{(r)}(r) + V^{(k)}(r) + V^{(s)}(r) + V^{(d)}(r)\ \ \ \ \ \ \ \ (4.4.1)
$$  

それぞれの項は、実空間で早期収束、逆空間で早期収束、自己エネルギー、双極子による補正項を表す。それぞれの中身は、  

$$  
V^{(r)}(r)= \frac{1}{2} \sum_n \sum^N_{i=1} \sum^N_{j=1}\frac{e^2}{4\pi\epsilon}\frac{Z_iZ_j}{|\boldsymbol{r_{ij}}+L\boldsymbol{n}|}\cdot erfc(\alpha |\boldsymbol{r_{ij}}+L\boldsymbol{n}|)\ \ \ \ \ \ \ \ (4.4.2)
$$  

$$  
V^{(k)}(r)= \frac{1}{2L^3}\sum_{\boldsymbol{k\not ={0}}} \frac{4\pi}{k^2}e^{-k^2/4\alpha^2}|\tilde{\rho}(\boldsymbol{k})|^2\ \ \ \ \ \ \ \ (4.4.3)
$$  

$$  
V^{(s)}(r)=-\frac{\alpha}{\sqrt{\pi}}\sum^N_{i=1} Z_i^2\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (4.4.4)
$$  

$$  
V^{(d)}(r)=\frac{2\pi}{(1+2\epsilon')L^3}\left( \sum^N_{i=1} Z_i\boldsymbol{r_i} \right)^2\ \ \ \ \ \ \ \ (4.4.5)
$$  

ここで、各記号の意味は、
  
|記号|単位|意味|
|:-:|:-:|:-:|
|L|Å|単位胞の長さ|
|$\alpha$|$L^{-1}$|Ewaldパラメータ|

また、  

$$  
erfc(r) = \frac{2}{\sqrt{\pi}} \int^{\infty}_r e^{-t^2}dt\ \ \ \ \ \ \ \ (4.4.6)
$$  

$$  
\boldsymbol{k} = \frac{2\pi \boldsymbol{m}}{L}\ \ \ \ \ \ \ \ (4.4.7)
$$  

$$  
\boldsymbol{n} = \left[
    \begin{matrix}
    n_x \\
    n_y \\
    n_z   
    \end{matrix}
\right]\ \ \ \ \ \ \ \ (4.4.8)
$$  

それぞれ、相補誤差関数、波数ベクトル、イメージセルのインデックスを表す整数ベクトルである。さらに、 $\tilde{\rho}(\boldsymbol{k})$ は電荷密度のフーリエ変換を表すので、  

$$  
\tilde{\rho}(\boldsymbol{k}) = \int_V d\boldsymbol{r} \rho(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{r}}=\sum^N_{j=1}Z_je^{-i\boldsymbol{k}\cdot\boldsymbol{r_j}} \ \ \ \ \ \ \ (4.4.9)
$$  

と表される。  
　では、エネルギーを分解したところで、偏微分してイオン $i$ にかかる力を求めよう。自己エネルギーは力に寄与しないので、  

$$  
\boldsymbol{F_i}=-\frac{\partial E^{(r)}}{\partial \boldsymbol{r_i}}-\frac{\partial E^{(k)}}{\partial \boldsymbol{r_i}}-\frac{\partial E^{(d)}}{\partial \boldsymbol{r_i}}\ \ \ \ \ \ \ (4.4.10)
$$  

さて、中身はどうなっているかというと、


$$  
-\frac{\partial E^{(r)}}{\partial \boldsymbol{r_i}}= \frac{e^2}{4\pi\epsilon} \sum_n \sum^N_{j=1}Z_iZ_j  \left(\frac{2\alpha}{\sqrt{\pi}} exp\left(-\alpha^2|\boldsymbol{r_{ij}}+L\boldsymbol{n}|^2\right)+\frac{erfc(\alpha |\boldsymbol{r_{ij}}+L\boldsymbol{n}|)}{|\boldsymbol{r_{ij}}+L\boldsymbol{n}|}\right)\frac{\boldsymbol{r_{ij}}+L\boldsymbol{n}}{|\boldsymbol{r_{ij}}+L\boldsymbol{n}|^2}\ \ \ \ \ \ \ \ (4.4.11)
$$  

$$  
-\frac{\partial E^{(k)}}{\partial \boldsymbol{r_i}}= \frac{1}{L^3}\sum_{j=1}^N Z_iZ_j \sum_{\boldsymbol{k\not ={0}}} \frac{4\pi\boldsymbol{k}}{k^2}e^{-k^2/4\alpha^2}sin(\boldsymbol{k\cdot r_{ij}})\ \ \ \ \ \ \ \ (4.4.12)
$$  

$$  
-\frac{\partial E^{(d)}}{\partial \boldsymbol{r_i}}=\frac{4\pi Z_i}{(1+2\epsilon')L^3} \sum^N_{j=1} Z_j\boldsymbol{r_j}\ \ \ \ \ \ \ \ (4.4.13)
$$  

このように分解し、エネルギーを早期に減衰させることで、エネルギーを計算している。

## 5. 構造最適化
MD計算をする際、初期構造が不安定なもの（原子間距離が近すぎてポテンシャルが大きいなど）のとき、温度が10万Kなどに上昇して計算が終了してしまうことがある。もちろん、毎回行わなければいけないわけではないが、行った方が良い。この構造最適化についても、いくつか方法があるが、このドキュメントでは __最急降下法__ と __共役勾配法__ について紹介する。

### 5.1 最急降下法

### 5.2 共役勾配法

## 6. 温度制御
まず、MD計算において温度がどのように定義されているかを説明する。系の温度Tは、  

$$  
T \equiv \frac{2}{3Nk_B} \sum_{i=1}^N \frac{1}{2}m_iv_i^2
$$  

とあらわされ、__運動学的温度__ と呼ばれる。熱力学で定義されている、__熱力学温度__ とは異なるので注意である。  
　温度がT一定の系に対して温度を制御する際は、系が温度Tの無限に大きい熱浴に浸かって熱平衡状態になっていると考えられる。このドキュメントでは、__Langevin法、能勢法、能勢-Hoover法__ の3つを紹介する。

### 6.1 Langevin法

### 6.2 能勢法

### 6.3 能勢-Hoover法

## 7. 圧力制御
まず、MD計算において圧力がどのように定義されているかを説明する。系の圧力Pは、  

$$  
P \equiv \frac{1}{3V}\left< \sum^N_{i=1} m\boldsymbol{v_i}^2 \right> + \frac{1}{3V}\left< \sum^N_{i=1} m\boldsymbol{r_i\cdot f_i}^{int} \right>\ \ \ \ \ \ (7.1)
$$  

$\boldsymbol{f}_i^{int}$ は原子間相互作用による力で、このPは、__内圧__ と呼ばれる。式(7.1)から、系の体積が運動すると考えれば圧力制御ができると考えた、__Andersen法__ を紹介したのち、立方体セルのみでなく任意の平行六面体にも適用可能にした __Parrinello-Rahman法__ を紹介して終わります。  

### 7.1 Andersen法

### 7.2 Parrinello-Rahman法