cmake ./
make

cp dtmc ./Examples/Higgs/
cp dtmc ./Examples/Micro/15m_rp/
cp dtmc ./Examples/Micro/email_cp/

cd ./Examples/Higgs/
tar xf Higgs.tar.xz
python plot_3_steps.py

cd ../Micro/15m_rp/
tar xf 15m_follower_s.net.tar.xz
python Micro_15m_rp.py

cd ../email_cp/
python Micro_email_cp.py

