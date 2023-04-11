# GK法のテスト用プロジェクト
main.f90のモジュール内では
サンプルコードはガウス関数に定数a,bを加えたものを積分しています。プログラム内では、a=b=0としています。ここを変えればいろんな関数に使えるはずです。
サブルーチン名の最後にwriteがついているものは被積分関数に対し、どのように刻んだかをdatファイルに書き出しながら積分するサブルーチンです。

gk法のソースコードはgk_test/gk.f90にまとめました。
これは、Notionの参考サイト内のシキノートさんのシンプルQUADPACKというタイトルのサイト内のコードをコピーし、module化したものです。

また、関数f(x)と積分区間と分点数を受け取り、中点法で数値積分する関数integral_midpointをmidpoint.f90にまとめました。

gk法のコードを全て理解するのは難しいです。引数が何かだけわかれば良いでしょう。
eps(in) : 推定誤差
ier(out) : 適応型積分の分割数
のはずです。
初めはmidpoint.f90をmain.f90から呼び出して使う方法を学んでイメージを掴んだら、その次にgk法を使うと良いでしょう。


## 注意点
ファイルの分割について
積分法のコード(gk.f90, midpoint.f90)と別の自作プログラムのコード(main.f90)（川村くんが使う場合はNSRのコードなど自分で作ったものになるはずです。）は、別のファイルに保存しています。通常、積分法のプログラムは、再利用ができるものです。すなわち、このファイルをプロジェクトのディレクトリにコピーすればNSRでもTMAでもT_c以下や絶対零度のコード、その他色んな数値積分が必要なプロジェクトで使い回せます。このようにプロジェクトによらないコードは機能ごとにファイルにまとめておくと便利です。また、便利なだけでなく、バグの混入を防ぐことが出来ます。
gk法のプログラムは2000行以上あるので、一つのファイルに入れると見づらいという理由もあります。
また、多くの便利な数値計算コードはネット上でライブラリとして共有されており、Githubなどからライブラリをダウンロードすることができます。
このような該当ファイルを使ってコードを組むことは多いです（Ex 高速フーリエ変換など）。


## makefileについて
このようにファイルを分割した時はifort gk.f90 midpoint.f90 main.f90 のように複数ファイルを依存関係に応じた順番でリンクしてコンパイルする必要があります。いちいちこれを毎回ターミナル上に打つのは面倒です。そこで登場するのがmakefileです。
gk_test/makefileの中身のように書くことで、上記コンパイル内容をあらかじめファイルにまとめておくことができます。
すぐに使えるように簡単に書いてありますが、編集する場合、makefileの書き方はネットなどで調べて見てください。

### 実行法
makefile内のあるディレクトリ内で、ターミナル上にmakeと打つと、実行ファイルtest(今回はtestで指定したが、いつもはデフォルトでa.outになります。)が生成されます。./testとターミナルに打てば実行されるはずです。


### シェルスクリプト
また、このようなコンパイルは、シェルスクリプトを用いて実行することも可能です。
シェルスクリプトは、ターミナル操作を自動化する際に用いるプログラミング言語のようなもので、より広い用途で使用されます。詳しくは調べてください。


## 多変数関数の一変数に着目し、それ以外の変数を固定した一変数関数を用意し、一変数関数を引数に持つサブルーチンに代入する方法
NSRプログラム内では、被積分関数は多変数関数のはずです（対相関関数Piなら、func_Pi(T, mu, q, nu, p)のように温度や化学ポテンシャル、重心エネルギー運動量などに依存する）。そして、この多変数関数の着目変数(Piなら相対運動量pに着目して積分しますよね)以外を固定して、着目変数の1変数関数として積分ルーチンに代入したいはずです。
これを行う方法が、mainプログラム内のモジュールmod_test_gkのサブルーチンintegrate_func...内に書いてあります。
このルーチンは多変数関数test_func(x, a, b)のa, bを固定し、xについて積分を行っています。
integrate_func内にcontains文の下に一変数関数であるラッパー関数wrapper(x)を用意します。
このラッパー関数の返り値retvalを多変数関数test_func(x, a, b)とします。
ここで、変数a, bについてはwrapper関数内で宣言しておらず、integrate_func内で宣言してことに注意してください。
実は、contains上の関数(integrate_func)で宣言した変数とcontains文以下の「関数内副プログラム」の変数は共有されるという性質があり、これを利用することで、多変数関数の変数a,bを関数内で設定した値で固定化し、着目変数xについての一変数関数をcontains以下に作ることができます。
そしてwrapper(x)はもはや一変数関数なので積分ルーチンdqag_k2()や、integral_midpoint()に代入することができます。
