# VASPに使用されている理論
参考文献は主に[VASPwiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)から。
## 1. 密度汎関数法
&emsp;VASPで用いられている第一原理計算の計算手法の一つである。

## 2. Smearing  
&emsp;INCARタグの __ISMEAR__ と __SIGMA__ に関して。  
　まず、電子の占有率（Fermi分布）などからエネルギーバンド計算や電子のSCF計算を行う際、エネルギーを離散的として計算すると、数値的な齟齬が生じてしまう。  
![SMEARING](https://github.com/MDGroup-WatanabeLab/image_for_mdpython/assets/138444525/a3ece222-b23b-444c-a534-fe1c577931b0)   
　そのため、上の図のオレンジ線色のように、エネルギーに幅を持たせ（SIGMA）、エネルギーを連続的であると近似することで、よりリアルに近い値をことができる。この幅を持たせることが「Smearing」である。

### 2.1. Gauss Smearing ( ISMEAR = 0 )

### 2.2. Fermi Smearing ( ISMEAR = -1 )

## 3. SCF 計算

## 4. 擬ポテンシャル法

## 5. 局所密度近似法