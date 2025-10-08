# QDays Prize

## メールアドレス

hk.quantum@icloud.com

## 経歴

日本のソフトウェア開発企業の会社員。  
量子コンピュータの仕組みを理解するため、自作の量子コンピュータシミュレータなどを作成。  
RSA暗号解読（2つの素数の合成数を素因数分解）では、自作の量子回路シミュレータを利用して、以下のビット数を解読する量子回路のシミュレーションに成功。

- ショアのアルゴリズムでは10bitを解読
- グローバーのアルゴリズムでは40bitを解読

## 今回解読した内容

量子回路をシミュレータによって実行できた最大ビット数は11bit。

```
--- Bit size 11 ---
Bit size: 11
Prime p: 1051
Curve order (#E): 1093
Subgroup order n: 1093
Cofactor h: 1
Generator point G: (471, 914)
Private key d: 756
Public key Q: (179, 86)
```

実機(ibm_torino)で実行できた最大ビット数は4bit。

```
--- Bit size 4 ---
Bit size: 4
Prime p: 13
Curve order (#E): 7
Subgroup order n: 7
Cofactor h: 1
Generator point G: (11, 5)
Private key d: 6
Public key Q: (11, 8)
```

## 実行環境

量子回路はqiskitで実装し、自作の高性能量子コンピュータシミュレータ `SparseStatevectorSimulator` をBackendとして実行。

IBMの量子コンピュータ「ibm_torino」では4bitのECC Curveを実行できたものの、ノイズのために結果はランダムとなり、期待通りの精度は得られなかった。  
なお、100ショットの実行時間は35秒程度であった。

## 実行手順

### 初期セットアップ

Python仮想環境を構築し、必要なライブラリをインストールする。

```
git clone https://github.com/hk-quantum/qday-prize.git
cd qday-prize
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

`.env`に、IBM Quamtum Platform で取得したAPIキーとCRNを設定する。

```
API_KEY=<Your API Key>
IBM_CRN=<Your CRN>
```

### 実行方法

起動引数に解読するビット数を指定する（例では11ビット）。

```
python src/main.py 11
```

`data/curves.json`から対象のビット数のECCパラメータを読み込み、量子回路構築およびシミュレータの実行ログを出力しつつ、最後に100ショットの測定結果に基づく秘密鍵の解読成功回数が出力される。

IBMの量子コンピュータ（実機）で実行する場合は以下のように起動する。

```
python src/ibm_main.py 3
```

### 実行結果

100ショットの測定結果から、解読した秘密鍵と解読が成功した回数を、ログの最後の行に出力している。

|log file|output result|backend|
|---|---|---|
|`logs/3_ibm.txt`|`Success: d=3 count=33`|ibm_torino|
|`logs/4_ibm.txt`|`Success: d=6 count=20`|ibm_torino|
|`logs/3.txt`|`Success: d=3 count=68`|SparseStatevectorSimulator|
|`logs/4.txt`|`Success: d=6 count=65`|SparseStatevectorSimulator|
|`logs/6.txt`|`Success: d=18 count=83`|SparseStatevectorSimulator|
|`logs/7.txt`|`Success: d=56 count=83`|SparseStatevectorSimulator|
|`logs/8.txt`|`Success: d=103 count=85`|SparseStatevectorSimulator|
|`logs/9.txt`|`Success: d=135 count=84`|SparseStatevectorSimulator|
|`logs/10.txt`|`Success: d=165 count=90`|SparseStatevectorSimulator|
|`logs/11.txt`|`Success: d=756 count=88`|SparseStatevectorSimulator|

