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

量子回路をシミュレータによって実行できた最大ビット数は12bit。

```
--- Bit size 12 ---
Bit size: 12
Prime p: 2089
Curve order (#E): 2143
Subgroup order n: 2143
Cofactor h: 1
Generator point G: (1417, 50)
Private key d: 1384
Public key Q: (1043, 1795)
```

実機(ibm_fez)で実行できた最大ビット数は5bit。

```
"bit_length": 5,
"prime": 23,
"a": 1,
"b": 11,
"generator_point": [7, 4],
"public_key": [9, 6]
```

## 実行環境

量子回路はqiskitで実装し、自作の高性能量子コンピュータシミュレータ `SparseStatevectorSimulator` をBackendとして実行。

IBMの量子コンピュータ「ibm_fez」では5bitのECC Curveを実行できたものの、量子回路が深くエラーのために結果はランダムとなり、期待通りの精度は得られなかった。  
なお、100ショットの実行時間は3bitでは3秒、4bitでは8秒、5bitでは30秒程度であった。

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

IBMの実機を利用する場合は、`.env`に IBM Quamtum Platform で取得したAPIキーとCRNを設定する。

```
API_KEY=<Your API Key>
IBM_CRN=<Your CRN>
```

### 実行方法

起動引数に解読するビット数とcompactまたはwideのアルゴリズムの種類を指定する。  
また、オプションとしてnumbaのJITコンパイラを利用する -single(単一スレッド) と -parallel(マルチスレッド)を用意している。

```
python src/main.py <num_bits> compact|wide [-single|-parallel]
```

`data/curves.json`から対象のビット数のECCパラメータを読み込み、量子回路構築およびシミュレータの実行ログを出力しつつ、最後に100ショットの測定結果に基づく秘密鍵の解読成功回数が出力される。

IBMの量子コンピュータ（実機）で実行する場合は以下のように起動する。

```
python src/ibm_main.py <num_bits> compact|wide
```

### 実行結果

100ショットの測定結果から、解読した秘密鍵と解読が成功した回数を、ログの最後の行に出力している。

以下は量子コンピュータの実機 `ibm_fez` で実行した結果を示す。

|log file|type|output result|
|---|---|---|
|`logs/3_ibm_compact.out`|compact|`Success: d=6 count=24`|
|`logs/3_ibm_wide.out`|wide|`Success: d=3 count=34`|
|`logs/4_ibm_compact.out`|compact|`Success: d=6 count=27`|
|`logs/4_ibm_wide.out`|wide|`Success: d=6 count=26`|
|`logs/5_ibm_compact.out`|compact|`Success: d=6 count=24`|

以下は自作の量子コンピュータシミュレータ `SparseStatevectorSimulator`で実行した [Version 4](https://github.com/hk-quantum/qday-prize/tree/v4.0.0)での結果である。

|log file|type|output result|
|---|---|---|
|`logs/3_compsact.out`|compact|`Success: d=3 count=68`|
|`logs/3_wide.out`|wide|`Success: d=3 count=75`|
|`logs/4_compact.out`|compact|`Success: d=6 count=79`|
|`logs/4_wide.out`|wide|`Success: d=6 count=75`|
|`logs/5_compact.out`|compact|`Success: d=6 count=80`|
|`logs/5_wide.out`|wide|`Success: d=6 count=75`|
|`logs/6_compact.out`|compact|`Success: d=18 count=82`|
|`logs/6_wide.out`|wide|`Success: d=18 count=76`|
|`logs/7_compact.out`|compact|`Success: d=56 count=85`|
|`logs/7_wide.out`|wide|`Success: d=56 count=85`|
|`logs/8_compact.out`|compact|`Success: d=103 count=87`|
|`logs/8_wide.out`|wide|`Success: d=103 count=88`|
|`logs/9_compact.out`|compact|`Success: d=135 count=87`|
|`logs/9_wide.out`|wide|`Success: d=135 count=85`|
|`logs/10_compact.out`|compact|`Success: d=165 count=88`|
|`logs/10_wide.out`|wide|`Success: d=165 count=95`|
|`logs/11_compact.out`|compact|`Success: d=756 count=92`|
|`logs/12_compact.out`|compact|`Success: d=1384 count=88`|

ビット数が多い区間では、シミュレータの性能限界により `compact` のみの実行結果を記載している。