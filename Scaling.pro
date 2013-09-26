PRO Scaling

  SET_PLOT, 'PS'
  DEVICE, /ENCAPSULATED    $
        , FILENAME         $
           = 'Scaling.eps' $
        , /COLOR           $
        , /INCHES          $
        , XSIZE=8.5        $
        , YSIZE=8.5        $
        , SCALE_FACTOR=1.0 $
        , XOFFSET=1.0      $
        , YOFFSET=1.0

  ; --- Strong Scaling (256^3 zones): 

  N_MPI_TASKS $
    = [16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0]

  ; End Time = 1

  WALL_TIME $
    = [108.93974900245667, $
       55.947239160537720, $
       28.285549163818359, $
       14.527857065200806, $
       7.6986460685729980, $
       3.9580910205841064, $
       1.8863770961761475, $
       0.99521517753601074]

  ; End Time = 3

  WALL_TIME $
    = [ $
       319.38515710830688, $
       164.23992800712585, $
       83.069401979446411, $
       42.278923988342285, $
       21.430653095245361, $
       11.023397922515869, $
       4.9565329551696777, $
       2.2640008926391602]

  PLOT, N_MPI_TASKS, WALL_TIME[0]/WALL_TIME $
      , LINE = 0, PSYM = - 6, SYMSIZE = 2   $
      , THICK = 6.0, CHARTHICK = 4.0        $
      , XR = [1.0E1, 3.0E3], XS = 1, /XLOG  $
      , XTITLE = 'MPI Tasks'                $
      , XTHICK = 6                          $
      , XCHARSIZE = 1.75                    $
      , XMARGIN = [14,4]                    $
      , YR = [0.8, 2.2E2], YS = 1, /YLOG    $
      , YTITLE = 'Speedup'                  $
      , YTHICK = 6                          $
      , YCHARSIZE = 1.75                    $
      , YMARGIN = [6,3]

  OPLOT, N_MPI_TASKS, N_MPI_TASKS/N_MPI_TASKS[0] $
       , LINE = 2, THICK = 6.0

  ; --- Weak Scaling (32^3 per MPI Task)

  N_MPI_TASKS $
    = [512, 4096, 13824]

  WALL_TIME $
    = [10.861274003982544, $
       18.460222005844116, $
       ]

  DEVICE, /CLOSE
  SET_PLOT, 'X'

END
