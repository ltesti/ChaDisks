import os
import time
import numpy as np

#  Function that reads the original ms table and writes output tables for the models and residual
#  NB. No checks done
def write_ascii_to_ms(origdata_ms,mod_txt,mod_ms,res_txt,res_ms,tb,dualpol=True,dores=True):
    print " "
    print "Create CASA ms with model ascii uvtable"
    print " "


    # Read ASCII file
    start=time.clock()

    f_in = open(mod_txt, "r")

    # Initialize arrays
    um = []; vm = []; Rem = []; Imm = []; Wem = [];

    # Go over all lines in ascii file and store data
    for line in f_in:

      # Store data in arrays
      variable = line.strip().split()

      um.append(float(variable[0]))
      vm.append(float(variable[1]))
      Rem.append(float(variable[2]))
      Imm.append(float(variable[3]))
      Wem.append(float(variable[4]))

    # Close ascii file
    f_in.close()

    if dores:
      f_in = open(res_txt, "r")

      # Initialize arrays
      ur = []; vr = []; Rer = []; Imr= []; Wer = [];

      # Go over all lines in ascii file and store data
      for line in f_in:

        # Store data in arrays
        variable = line.strip().split()

        ur.append(float(variable[0]))
        vr.append(float(variable[1]))
        Rer.append(float(variable[2]))
        Imr.append(float(variable[3]))
        Wer.append(float(variable[4]))

      # Close ascii file
      f_in.close()

    print "Finished reading model visibilities..."
    end = time.clock()
    print "It took", end-start, "seconds"

    # Read ms table and replace data column
    start=time.clock()


    # Create new tables
    syscommand = 'rm -rf '+ mod_ms
    os.system(syscommand)
    syscommand = 'cp -r '+ origdata_ms + ' ' + mod_ms
    os.system(syscommand)

    if dores:
      syscommand = 'rm -rf '+ res_ms
      os.system(syscommand)
      syscommand = 'cp -r '+ origdata_ms + ' ' + res_ms
      os.system(syscommand)

    # Open ms table
    tb.open(mod_ms)

    # Check if CORRECTED_DATA column is present:
    # (MODEL_DATA column will be modified when cleaning)
    all_columns = tb.colnames()
    correct     = False
    if 'CORRECTED_DATA' in all_columns:
      correct_data_m = tb.getcol("CORRECTED_DATA")
      correct = True

    # Get column UVW, DATA, WEIGHT and DATA_DESC_ID (spw info!)
    data_m   = tb.getcol("DATA")
    uvw_m    = tb.getcol("UVW")
    weight_m = tb.getcol("WEIGHT")
    spw_m    = tb.getcol("DATA_DESC_ID")
    tb.close()

    if dores:
      tb.open(res_ms)

      # Check if CORRECTED_DATA column is present:
      # (MODEL_DATA column will be modified when cleaning)
      all_columns = tb.colnames()
      correct     = False
      if 'CORRECTED_DATA' in all_columns:
        correct_data_r = tb.getcol("CORRECTED_DATA")
        correct = True

      # Get column UVW, DATA, WEIGHT and DATA_DESC_ID (spw info!)
      data_r   = tb.getcol("DATA")
      uvw_r    = tb.getcol("UVW")
      weight_r = tb.getcol("WEIGHT")
      spw_r    = tb.getcol("DATA_DESC_ID")
      tb.close()


    # Go over msdata array and input new model values into
    # the XX (data[0,0,:]) and YY (data[1,0,:]) components:
    print "Size of current ms file: ", len(spw_m)
    for i in range(len(spw_m)):

      data_m[0,0,i] = complex(Rem[i],Imm[i])
      weight_m[0,i] = Wem[i]
      if dualpol:
        data_m[1,0,i] = complex(Rem[i],Imm[i])
        weight_m[1,i] = Wem[i]
      if correct is True:
        correct_data_m[0,0,i] = complex(Rem[i],Imm[i])
        if dualpol:
          correct_data_m[1,0,i] = complex(Rem[i],Imm[i])
      if dores:
        data_r[0,0,i] = complex(Rer[i],Imr[i])
        weight_r[0,i] = Wer[i]
        if dualpol:
          data_r[1,0,i] = complex(Rer[i],Imr[i])
          weight_r[1,i] = Wer[i]
        if correct is True:
          correct_data_r[0,0,i] = complex(Rer[i],Imr[i])
          if dualpol:
            correct_data_r[1,0,i] = complex(Rer[i],Imr[i])

    # Put back modified weight column on original ms

    # Put back modified column
    tb.open(mod_ms,nomodify=False)
    tb.putcol("DATA",data_m)
    tb.putcol("WEIGHT",weight_m)
    if correct is True:
      tb.putcol("CORRECTED_DATA",correct_data_m)
    tb.flush()
    # Close table tools
    tb.close()

    if dores:
      tb.open(res_ms,nomodify=False)
      tb.putcol("DATA",data_r)
      tb.putcol("WEIGHT",weight_r)
      if correct is True:
        tb.putcol("CORRECTED_DATA",correct_data_r)
      tb.flush()

    # Close table tools
    tb.close()

    print "Finished with model     ms: ", mod_ms
    if dores:
      print "     and with residuals ms: ", res_ms
    end = time.clock()
    print "It took", end-start, "seconds"
    print ""


