# QDays Prize

## メールアドレス

hk.quantum@icloud.com

## 経歴

日本のソフトウェア開発企業の会社員。  
量子コンピュータの仕組みを理解するため、自作の量子コンピュータシミュレータなどを作成。  
RSA暗号解読では、自作の量子回路シミュレータを利用して、2つの素数を掛け算した合成数を素因数分解する量子回路を過去に作成。

- ショアのアルゴリズムでは10bitを解読
- グローバーのアルゴリズムでは40bitを解読

## 今回解読した内容

以下の11bitのECCを解読。

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

## 実行環境

量子回路はqiskitで実装し、自作の高性能量子コンピュータシミュレータ `SparseStatevectorSimulator` をBackendとして実行。

IBMの量子コンピュータ「ibm_torino」で実行を試みたものの、4bitであっても利用可能な cz, sx の量子ゲート上限を超えていたため実行できなかった。

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

### 実行方法

起動引数に解読するビット数を指定する（例では11ビット）。

```
python src/main.py 11
```

`data/curves.json`から対象のビット数のECCパラメータを読み込み、量子回路構築およびシミュレータの実行ログを出力しつつ、最後に100ショットの測定結果に基づく秘密鍵の解読成功回数が出力される。

### 実行結果

- `logs/11.log`: 11bitのシミュレーション時のログ
