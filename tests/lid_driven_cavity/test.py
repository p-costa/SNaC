#!/usr/bin/env python
def test_ldc():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    data_1,xp_1,yp_1,zp_1,xu_1,yv_1,zw_1 = read_single_field_binary("vey_fld_0001500_b_001.bin",1,np.array([1,1,1]))
    data_2,xp_2,yp_2,zp_2,xu_2,yv_2,zw_2 = read_single_field_binary("vey_fld_0001500_b_002.bin",2,np.array([1,1,1]))
    data_3,xp_3,yp_3,zp_3,xu_3,yv_3,zw_3 = read_single_field_binary("vey_fld_0001500_b_003.bin",3,np.array([1,1,1]))
    data_4,xp_4,yp_4,zp_4,xu_4,yv_4,zw_4 = read_single_field_binary("vey_fld_0001500_b_004.bin",4,np.array([1,1,1]))
    data_1 = data_1[0,:,:]
    data_2 = data_2[0,:,:]
    data_3 = data_3[0,:,:]
    data_4 = data_4[0,:,:]
    data = np.block([[data_1,data_3],[data_2,data_4]])
    data_ref = np.loadtxt("data_ldc_re1000.txt")
    islice = int(np.size(data[0,:])/2)
    np.testing.assert_allclose(data[islice,:], data_ref[:,1], rtol=1e-7, atol=0)
if __name__ == "__main__":
    test_ldc()
    print("Passed!")
