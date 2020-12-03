
if(exist("n")==0 || n<0) n = n0 # ループ変数の初期化

filename = sprintf("./output_Ez/%05d.dat", n) # n番目のデータファイルの名前の生成
time = sprintf("%.3f [ps] ", n*2.072)

unset label
set label time at graph 0.8,0.9

splot filename u 1:2:3 with pm3d

pause 0.2
#-------------------------------------------------------------------------------
# ループ処理
#-------------------------------------------------------------------------------
n = n + dn            # ループ変数の増加
if ( n < n1 ) reread  # ループの評価
undefine n            # ループ変数の削除
