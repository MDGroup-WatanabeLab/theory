# VASPに使用されている理論
参考文献は主に[VASPwiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)から。また、[解説スライド](https://www.slideshare.net/dc1394/ss-26378208)も参考にした。
<div style="text-align: right;">文責：内藤</div>

## 0. 必要になりそうな前提知識
　この後の話を理解するのに必要そうな式や用語をまとめておく。
### 0.1. Schrödinger方程式
　1926年にErwin R.J.A. Schrödinger が発表？した。式である。単一粒子かつ定常状態の場合は、  
$H^{\wedge}$


### 0.2. Dirac方程式

## 1. 密度汎関数法
&emsp;VASPで用いられている第一原理計算の計算手法の一つである。その概要について述べる。  
　DFT計算では、電子軌道内を運動する電子を、電子の密度分布で表現する。この電子密度と個々の電子間に働く引力・斥力相互作用を考慮することで、系内に存在する全ての電子間の相互作用を評価できるため、各原子・分子の安定構造を計算することが可能である。（[参考](https://rdreview.jaea.go.jp/review_jp/kaisetsu/723.html#:~:text=DFT%E8%A8%88%E7%AE%97%E3%81%A7%E3%81%AF%E3%80%81%E9%9B%BB%E5%AD%90%E8%BB%8C%E9%81%93,%E3%81%93%E3%81%A8%E3%81%8C%E5%8F%AF%E8%83%BD%E3%81%A7%E3%81%82%E3%82%8B%E3%80%82)）  
　すなわち、MD計算では、各原子間に働く相互作用を原子間ポテンシャルからエネルギー・力を求めることで軌跡を計算していたが、DFT計算では原子ではなく __電子__ について相互作用を考えていくということである。
　

## 2. Smearing  
&emsp;INCARタグの __ISMEAR__ と __SIGMA__ に関して。  
　まず、電子の占有率（Fermi分布）などからエネルギーバンド計算や電子のSCF計算を行う際、エネルギーを離散的として計算すると、数値的な齟齬が生じてしまう。  
![SMEARING](https://github.com/MDGroup-WatanabeLab/image_for_mdpython/assets/138444525/a3ece222-b23b-444c-a534-fe1c577931b0)   
　そのため、上の図のオレンジ線色のように、エネルギーに幅を持たせ（SIGMA）、エネルギーを連続的であると近似することで、よりリアルに近い値をことができる。この幅を持たせることが「Smearing」である。

### 2.1. Gauss Smearing ( ISMEAR = 0 )

### 2.2. Fermi Smearing ( ISMEAR = -1 )

## 3. SCF 計算

### 3.1. The iterative subspace (RMM-DIIS)

### 3.2. blocked Davidson algorithm

## 4. 擬ポテンシャル法

## 5. 局所密度近似法