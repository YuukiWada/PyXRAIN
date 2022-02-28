#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import binascii
import math
import time

class xrain:
    def __init__(self,input_file,switch=True):
        self.par = self.read_par(input_file)
        self.ppi = self.read_ppi(input_file)
        self.sect = self.read_sect(input_file)
        self.lat = round(self.par[11][0]+(self.par[11][1]/60.0)+(self.par[11][2]/3600.0),4)
        self.lon = round(self.par[12][0]+(self.par[12][1]/60.0)+(self.par[12][2]/3600.0),4)
        self.alt = self.par[13]
        self.az_start = self.sect[0]
        self.az_end = self.sect[1]
        self.el_start = self.sect[2]
        self.el_end = self.sect[3]
        self.range_num = self.par[22]
        self.range_min = self.par[19]
        self.range_max = self.par[20]
        self.range_step = self.par[21]
        self.az_num = self.par[23]
        self.date = self.par[2][0]+"-"+self.par[2][1]+"-"+self.par[2][2]
        self.time = self.par[3][0]+":"+self.par[3][1]+":00"
        self.obs_start = self.par[17][0]+":"+self.par[17][1]+":"+self.par[17][2]
        self.obs_end = self.par[18][0]+":"+self.par[18][1]+":"+self.par[18][2]
        self.elv = self.par[9]
        self.x = self.axis_x()
        self.y = self.axis_y()
        self.z = self.axis_z()
        self.rng = self.axis_range()
        if switch==True:
            self.credit()

    def credit(self):
        print("")
        print("  ------------------------------------------------------")
        print("     国土交通省 XRAIN バイナリ解析 Python ライブラリ")
        print("     Version 1.0 (2022年2月28日)")
        print("     作成: 和田有希 (大阪大学大学院工学研究科)")
        print("     https://github.com/YuukiWada/PyXRAIN")
        print("")
        print("     本ソフトウェアを用いて生じた不利益・損害について")
        print("     作成者は一切の責任を負いません。")
        print("  ------------------------------------------------------")
        print("")

    def read_par(self,input_file):
        f = open(input_file, "rb")
        data = f.read()
        par = list()
        par.append(binascii.hexlify(data[4:6]))                       # 0: observation site
        par.append(binascii.hexlify(data[7]))                         # 1: data type
        par.append([data[8:12],data[13:15],data[16:18]])              # 2: observation date
        par.append([data[19:21],data[22:24]])                         # 3: observation time
        if binascii.hexlify(data[33])=="01":                          # 4: status code
            par.append(0)
        else:
            par.append(1)
        par.append(float(binascii.hexlify(data[40:42]))/10.0)         # 5: rotation speed
        if binascii.hexlify(data[42:44])=="0001":                     # 6: identification of PPI/CAPPI
            par.append("CAPPI")
        else:
            par.append("PPI")
        par.append(int(binascii.hexlify(data[44:46]),16))             # 7: elevation step
        par.append(int(binascii.hexlify(data[46:48]),16))             # 8: elevation step num
        if int(binascii.hexlify(data[48:50]),16)>32767:               # 9: elevation angle
            par.append((int(binascii.hexlify(data[48:50]),16)-2**16)/100.0)
        else:
            par.append(int(binascii.hexlify(data[48:50]),16)/100.0)
        if binascii.hexlify(data[52:56])=="00000004":                 # 10: site status
            par.append(0)
        else:
            par.append(1)
        par.append([int(binascii.hexlify(data[62:64]),16),int(binascii.hexlify(data[64:66]),16),int(binascii.hexlify(data[66:68]),16)]) # 11: latitude
        par.append([int(binascii.hexlify(data[68:70]),16),int(binascii.hexlify(data[70:72]),16),int(binascii.hexlify(data[72:74]),16)]) # 12: longitude
        par.append(int(binascii.hexlify(data[74:78]),16)/100.0)       # 13: altitude (m)
        par.append(int(binascii.hexlify(data[88:90]),16)/100.0)       # 14: horizontal power (kW)
        par.append(int(binascii.hexlify(data[102:104]),16)/100.0)     # 15: vertical power (kW)
        par.append(int(binascii.hexlify(data[110:112]),16))           # 16: transmission frequency (MHz)
        par.append([data[128:130],data[131:133],data[134:136]])       # 17: observation start time
        par.append([data[136:138],data[139:141],data[142:144]])       # 18: observation start time
        par.append(int(binascii.hexlify(data[144:148]),16)/100.0)     # 19: range start (m)
        par.append(int(binascii.hexlify(data[148:152]),16)/100.0)     # 20: range end (m)
        par.append(int(binascii.hexlify(data[152:156]),16)/100.0)     # 21: range resolution (m)
        par.append(int(binascii.hexlify(data[156:160]),16))           # 22: range num
        par.append(int(binascii.hexlify(data[160:162]),16))           # 23: azimuth num
        f.close()
        return par

    def read_ppi(self,input_file):
        par = self.par
        f = open(input_file, "rb")
        data = f.read()
        num_rng = par[22]
        num_az = par[23]
        value = np.zeros((num_az, num_rng))
        for i in range(num_az):
            for j in range(num_rng):
                index = 512+i*(16+2*num_rng)+(16+2*j)
                value[i,j] = int(binascii.hexlify(data[index:index+2]),16)
        data_type = par[1]
        type_num = ["05", "06", "07", "08", "09", "0E", "11", "12", "15", "19", "21", "25", "31", "35"]
        a = [90.0, 95.0, 100.0, 105.0, 1.0, 80.0, 85.0, 1.0, 1.0, 1.0, 1.0, 1.0, 360.0, 1.0]
        b = [0.0, 0.0, 0.0, 0.0, 32768.0, 0.0, 0.0, 32768.0, 32768.0, 1.0, 32768.0, 1.0, 1.0, 32768.0]
        c = [16384, 16384, 16384, 16384, 100.0, 16384, 16384, 100.0, 100.0, 100.0, 100.0, 65533.0, 65534.0, 100.0]
        type_index = type_num.index(data_type)
        value = a[type_index]*(value-b[type_index])/c[type_index]
        f.close()
        return value

    def read_sect(self,input_file):
        num_rng = self.par[22]
        num_az = self.par[23]
        sect_info=list()
        f = open(input_file, "rb")
        data = f.read()
        for i in range(num_az):
            index = 512+i*(16+2*num_rng)
            sect_info.append([int(binascii.hexlify(data[index:index+2]),16)/100.0,int(binascii.hexlify(data[index+2:index+4]),16)/100.0,\
                              int(binascii.hexlify(data[index+4:index+6]),16)/100.0,int(binascii.hexlify(data[index+6:index+8]),16)/100.0])
        f.close()
        return sect_info

    def axis_x(self):
        n_rng = self.range_num
        n_az = self.az_num
        min_rng = self.range_min
        step_rng = self.range_step
        step_az = 360.0/n_az
        rng = np.tile(min_rng+np.arange(n_rng+1)*step_rng/1000.0, (n_az+1,1))
        az = np.tile(np.arange(n_az+1)*step_az, (n_rng+1,1)).transpose()
        x = rng*np.sin(az*math.pi/180)
        return x

    def axis_y(self):
        n_rng = self.range_num
        n_az = self.az_num
        min_rng = self.range_min
        step_rng = self.range_step
        step_az = 360.0/n_az
        rng = np.tile(min_rng+np.arange(n_rng+1)*step_rng/1000.0, (n_az+1,1))
        az = np.tile(np.arange(n_az+1)*step_az, (n_rng+1,1)).transpose()
        y = rng*np.cos(az*math.pi/180)
        return y

    def axis_z(self):
        n_rng = self.range_num
        min_rng = self.range_min
        step_rng = self.range_step
        rng = min_rng+(np.arange(n_rng)+0.5)*step_rng/1000.0
        elv = self.elv
        alt = self.alt
        z = rng*np.tan(elv*math.pi/180.0)+alt/1000.0
        return rng
    
    def axis_range(self):
        n_rng = self.range_num
        min_rng = self.range_min
        step_rng = self.range_step
        rng = min_rng+(np.arange(n_rng)+0.5)*step_rng/1000.0
        return rng
    
    def site(self):
        site_id = ["8105", "8106", "8107", "8108", "8109", "8205", "8206", "8207", "8208", "8209", "820a", "820b",\
                   "8305", "8306", "8405", "8406", "8407", "8408", "8409", "840a", "8505", "8506", "8507", "8508",\
                   "8605", "8606", "8607", "8608", "8609", "860a", "860b", "8705", "8706", "8707", "8708",\
                   "8805", "8806", "8807", "8808", "8900", "8a00"]

        site_name = ["関東/関東", "関東/新横浜", "関東/氏家", "関東/八丈島", "関東/船橋", "九州/九千地", "九州/菅岳", "九州/古月山", "九州/風師山", "九州/桜島", "九州/山鹿", "九州/宇城",\
                     "北海道/北広島", "北海道/石狩", "東北/一関", "東北/一迫", "東北/涌谷", "東北/岩沼", "東北/伊達", "東北/田村", "北陸/水橋", "北陸/能美", "東北/京ヶ瀬", "東北/中の口",\
                     "中部/尾西", "中部/安城", "中部/鈴鹿", "中部/静岡北", "中部/香貫山", "中部/富士宮", "中部/浜松", "近畿/六甲", "近畿/葛城", "近畿/鷲峰山", "近畿/田口",\
                     "中国/熊山", "中国/常山", "中国/野貝原", "中国/牛尾山", "四国", "沖縄"]
        index = site_id.index(self.par[0])
        if index==None:
            return "未定義"
        else:
            return site_name[index]

    def obs_type(self):
        obs_id = ["05", "06", "07", "08","09", "0E", "11", "12", "15", "19", "21", "25", "31", "35"]
        type_name = ["受信電力強度 (dB)", "受信電力強度 (dB)", "受信電力強度 (dB)", "受信電力強度 (dB)", "受信電力強度 (dB)", "受信電力強度 (dB)", "受信電力強度 (dB)",\
                     "レーダー反射因子 (dBZ)", "風速 (m/s)", "分散 (m/s)", "反射因子差 (dB)", "偏波間相関係数", "偏波間位相差 (deg)", "伝播位相差変化率 (deg/km)"]
        index = obs_id.index(self.par[1])
        return type_name[index]

    def parameter(self):
        par = self.par
        status = ["正常","異常"]
        message = list()
        lat = round(par[11][0]+(par[11][1]/60.0)+(par[11][2]/3600.0),4)
        lon = round(par[12][0]+(par[12][1]/60.0)+(par[12][2]/3600.0),4)
        message.append("   レーダーサイト     : "+str(self.site()))
        message.append("   データ種別         : "+str(self.obs_type()))
        message.append("   データ配信日時     : "+str(par[2][0])+"年"+str(par[2][1])+"月"+str(par[2][2])+"日 "+str(par[3][0])+"時"+str(par[3][1])+"分 (日本標準時)")
        message.append("   XRAINステータス    : "+status[par[4]])
        message.append("   サイトステータス   : "+status[par[10]])
        message.append("   レーダー回転速度   : "+str(par[5])+" rpm")
        message.append("   運用モード         : "+par[6])
        message.append("   CAPPI仰角数        : "+str(par[7]))
        message.append("   CAPPI仰角番号      : "+str(par[8]))
        message.append("   CAPPI仰角          : "+str(par[7])+"度")
        message.append("   レーダーサイト緯度 : 北緯 "+str(lat)+"度")
        message.append("   レーダーサイト経度 : 東経"+str(lon)+"度")
        message.append("   レーダーサイト高度 : "+str(par[13])+" m")
        message.append("   水平偏波送信電力   : "+str(par[14])+" kW")
        message.append("   垂直偏波送信電力   : "+str(par[15])+" kW")
        message.append("   送信周波数         : "+str(par[16]/1000)+" GHz")
        message.append("   観測開始時刻       : "+str(par[17][0])+"時"+str(par[17][1])+"分"+str(par[17][2])+"秒 (日本標準時)")
        message.append("   観測終了時刻       : "+str(par[18][0])+"時"+str(par[18][1])+"分"+str(par[18][2])+"秒 (日本標準時)")
        message.append("   視線方向最小距離   : "+str(par[19])+" m")
        message.append("   視線方向最大距離   : "+str(par[20])+" m")
        message.append("   視線方向ステップ   : "+str(par[21])+" m")
        message.append("   視線方向観測数     : "+str(par[22]))
        message.append("   方位角観測数       : "+str(par[23]))
        for string in message:
            print(string)


