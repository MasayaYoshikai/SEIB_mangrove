�������ӁF�@���p������readme�������̂̓V���h�C�̂ŁA���̓��{���readme�ɂ��ẮASEIB-v2.61�ȍ~�̍X�V����߂܂����B
���̃t�@�C���́B�ꉞ�A�܂��Q�l�ɂȂ镔�������邩�Ǝv���Y�t���Ă��邾���ł��̂ŁA�R�[�h�̎g�p�ɍۂ��Ă�reade.txt���Q�Ƃ���ĉ������B

��������������������������������������������������������������������
SEIB-DGVM ��舵�������� (���{��ŁA2010�N11��30������)

�����i�i���É���w���w�����ȁj
http://seib-dgvm.com/hsato/

�@�@�@--- �ڎ� ---
�E�g�p��������
�E�t�@�C���\��
�E���s�菇�i���n����̎��s�j
�E���s�菇�i���X�^�[�g���s�j
�E�ݒ�t�@�C��
�ESEIB-Viewer
�Emain.f90���̃R�[�h�̗���
�E��̓v���O����
�E�G�f�B�^�ɂ���
�E���s�R�[�h�ɂ���
�E���̑��̒���
��������������������������������������������������������������������


���������@�g�p�������� (ver 1.0)�@��������
(1) SEIB-DGVM�̃R�[�h�́A�����p�r�Ɍ���l�E�@�l��킸�䗘�p���������܂��B�c����s���v��ւ̗��p����]�����ꍇ�ɂ́A�\�ߍ����i�̋��𓾂Ă��������B
(2) �R�[�h�̐��m����o�͂̐��x�Ɋւ��Ă̕ۏ؂͂��������˂܂��̂ŁA���p�ɂ���āA�����Ȃ��肪�������Ă��A�����ł͂��������ӔC�������܂���B
(3) �R�[�h�̍Ĕz�z�͎��R�ł����A���̍ہA�t�@�C���\������e�͕ύX���Ȃ��ŉ������B�܂��A���ς����R�[�h�̍Ĕz�z����]�����ꍇ�ɂ́A�\�ߍ����i�̋��𓾂Ă��������B
(4) ���̃R�[�h�𗘗p�����������ʂ��o�ł����ۂɂ́A�K�؂ȕ��������p���Ă��������B�Ȃ��A���̃R�[�h���g�p�����o�ŕ��ɂ��܂��ẮA�ꕔ������肦��΍K���ɑ����܂��B


���������@�t�@�C���\���@��������
���������́@���g�b�v�t�H���_
readme.txt  
history.txt  �X�V����

readme_j.txt   ���{��ɂ��R�[�h���p�@�i���̕��́j
readme_e.txt   �p��ɂ��R�[�h���p�@
history.txt    �X�V����
Model_description_??????_.pdf   �ŐV�Ń��f���̋L�q

���\�[�X�R�[�h�@���t�H���_"Code"
start_point.f90     �p�����[�^�[�ƋC�ۓy��f�[�^�̓ǂݏo��
modules.f90         �S���[�`���ŋ��L����萔�E�֐��̐錾
main.f90            ���C�����[�`��
initialize.f90      ������
physics.f90         �����ߒ��i���ˎ��x�A�����x�Ȃǁj
metabolic.f90       �����ߒ��i�������A�ċz�A�����Ȃǁj
spacial_calc.f90    �O������Ԍv�Z�i�̖��̌����Ȃǁj
population_regu.f90 �A�����ԉߒ��i���S�A�蒅�A�����Ȃǁj
output.f90          �o�̓t�@�C���̍쐬
etc.f90             ��̃J�e�S���[�ɓ���Ȃ����T�u�v���O����

���p�����[�^�[��`�t�@�C���@���t�H���_"Code"
parameter.txt

���T���v���C�ۃf�[�^�@���t�H���_"Code"
climate.txt     ���f�����͗p�̃T���v���C�ۃf�[�^�B�J���}��؂�̃e�L�X�g�`���BSEIB-DGVM�E�F�u�T�C�g��̋C�ۃf�[�^�W�F�l���[�^�[�ŁA�k��2�x58���A���o102�x18���A����1901�`2005�N���w�肵�č쐬�B�Ȃ��A���̈ʒu�̓}���[�V�A�̃p�\�����тɑ�������B

���y�n���f�[�^�@���t�H���_"Code"
land_prop.txt     ���f�����͗p�̑S���P�x�O���b�h���b�V���̓y�n���f�[�^�B�J���}��؂�̃e�L�X�g�`���B�k��90�x���o180�x���N�_�ɓ������Ƀf�[�^�����ׂ��A���o180�x�ɒB�����玟�͖k��89�x180�x���f�[�^���n�܂�BGlobal Soil Wetness Project 2�ihttp://grads.iges.org/gswp/�j�̓��͗p�f�[�^���_�E�����[�h���āA���`�������́B���p���@��f�[�^�̗��p�K���ɂ��ẮA����WEB�y�[�W�ɋL�q����������Ȃ������B�e�s�̃f�[�^�̕��я��͎��̒ʂ�FLand Mask (�C:0, ��:1)�A���ύ��x (m)�A�y��A���x�h�A�ő�ې���(saturation point)�A�ޏ�e����(field capacity)�A��z����(matrix potential)�A�ނ�_(wilting point)�B

����{��̓v���O�����@���t�H���_"Tools"
analysis1.f90        �X�ѓ��Ԃ̊�{�l���o�͂���v���O����
Viewer_snapshot.pov  ���z�ѕ��\���p��POV-RAY�v���O����(�Î~��)
Viewer_animation.pov  ���z�ѕ��\���p��POV-RAY�v���O����(�A����)

��SEIB-Viewer�@���t�H���_"SEIB-Viewer"
���̃t�H���_�ɂ́A�v���O�����{�̂ƁA���̎��s�ɕK�v�ȃt�@�C�����܂ނQ�̃T�u�t�H���_���[�߂��Ă��܂��B


���������@���s�菇�i���n����X�^�[�g����ꍇ�j�@��������
(1) ���W���w�肷��B star_point.f90����
  LAT     =  2.97 !north:+, south:- (decimalized)
  LON     = 102.3 ! east:+,  west:- (decimalized)

�Ƃ����ӏ����A�V�~�����[�g������ꏊ�̈ܓx�o�x�ɏ���������B�ܓx��+90�`-90�܂ł̐��l�����͉\�Ŗk�܂�+��܂�-�Ŏw��A�o�x��+180�`-180�����͉\�œ��o��+���o��-�Ŏw��B�Ȃ��A�����_�ȉ���"��"�ł͂Ȃ��ď\�i�@�ɂċL�����邱�ƁB

(2) �K�v�ł����modules.f90���̃p�����[�^�[��ς���B�����ɂ́A�Ⴆ�Έȉ��̂悤�ȃp�����[�^�[���܂܂�Ă���B
Max_no            �V�~�����[�g����ؖ{�̂̍ő吔
Max_hgt           �V�~�����[�g���鉼�z�ѕ��̍������C���[���B�����z�v�Z���s�킹��ۂɎQ�Ƃ���B�P���C���[�̌��݂́A�ȉ��Ő�������p�����[�^�[STEP�ɂ���`�����B�����]�T�������āi�ő�����{���j�ݒ肷�邱�ƁB
Dived             �ؖ{�蒅���b�V���ׂ̍����B�Ⴆ��10�Ǝw�肷��ƁA���z�ѕ��̗я���10�~10=100�ɋ敪���A�e�敪�̒��S�ɖؖ{���蒅�ł��鎖�ɂȂ�B���߂�Max_loc�Ɠ������x�̒l�ɂ��Ȃ��ƁA�X�т̃V�~�����[�V�����Ƃ��Ă͕s���R�B
DivedG             �я��ɑ��݂��鑐�{�Z���̉𑜓x�B�Ⴆ��10�Ǝw�肷��ƁA���z�ѕ��̗я���10�~10=100�̑��{�Z�������݂��邱�ƂɂȂ�B�e���{�Z���͓Ɨ��ɁA���ꂼ��̏ꏊ�̌������ɂ����鐬���̃V�~�����[�V�������ɍs���B�A���A���ȊO�̊������i�y��ܐ����Ȃǁj�͓������z�ѕ��̑S�Ă̑��{�Z���ŋ��ʂł���B

(3) ���s�t�@�C�����쐬����Bstart_point.f90����include���ő��̃v���O�����t�@�C�������w�肳��Ă���̂ŁAstart_point.f90�݂̂��R���p�C������悤�R�}���h��łĂ΁A�S�R�[�h�������I�ɕϊ������B��̓I�ɂ́A���s�t�@�C������go.out�Ƃ���̂Ȃ�A�R�}���h���C�����ȉ��̂悤�ɓ��͂���B
�@�@�@f90 start_point.f90 -o go.out
���̗�́A�R���p�C���R�}���h���uf90�v�̏ꍇ�BIntel Fortran�ł́uifort�v�Ag95�R���p�C���ł́ug95�v�Ƃ��邱�ƁB

(4) �p�����[�^�[�t�@�C��parameter.txt�̐ݒ�B�d�v�ȓ_�͎��̂Q�F(1)�C�ۃf�[�^�Ɠy�n�f�[�^�̒u����Ă���f�B���N�g����FullPath�Ŏw��B(2)Simulation_year�ɁA�V�~�����[�V����������N�����w�肷��B���s���̌v�Z���x���悭������Ȃ������́A20�N���x���w�肵�ėl�q�����邱�ƁB�Ȃ��Aparameter.txt�͎��s�t�@�C���Ɠ���f�B���N�g���ɒu�����ƁB

(5) ���s�t�@�C�����i��̗�Ȃ�go.out�AUnix��ł�./go.out�j����͂��āA�v���O���������s������B���s�ア�����̏o�̓t�@�C���ioutput_*.txt�j���쐬�����B��������āA�p�����[�^�[��R�[�h��ς��čēx���s������A���邢�͉�̓f�[�^�Ƃ��ė��p����A�Ƃ�������B�����o�̓f�[�^�̉�͂ɂ́A�Y�t������̓v���O�������֗��ł���B��̓v���O�����̎g�p�@�Ɋւ��Ă͌�q�B

(6) ��L�̎菇�ɂăv���O���������s�����̂��m�F������A�C�ۃf�[�^'climate.txt'���V�~�����[�g���������ꏊ�̃f�[�^�ƍ����ւ��邱�ƁB�C�ۃf�[�^�̓f�[�^��SEIB-DGVM�E�F�u�T�C�g�ɂĐ����\�B


���������@���s�菇�i�O��̑�������X�^�[�g����ꍇ�j�@��������
(1) ���X�^�[�g�f�[�^�̗p��
parameter.txt����Flag_spinup_write���A�u.true.�v�Ɛݒ肵�������ŁA��Ŏ��������@�ɏ]���ăV�~�����[�V���������s����B���s��A���X�^�[�g�p�f�[�^�ł���Aspinup_out.txt���o�͂����B
(2)���X�^�[�g
����ꂽspinup_out.txt���Aspinup_in.txt�ƃt�@�C������ς���B���̃t�@�C���́A���s�t�@�C���Ɠ����t�H���_�ɂ����K�v������Bparameter.txt����Flag_spinup_read���A�u.true.�v�Ɛݒ肵�������ŁA��Ŏ��������@�ɏ]���ăV�~�����[�V���������s����B���X�^�[�g��́Aparameter.txt����Simulation_year�Ŏw�肳���N�����A�V���ɃV�~�����[�g�����B


���������@�ݒ�t�@�C���@��������
�ݒ�t�@�C��(parameter.txt)�́A�R�[�h���s���ɓǂݍ��܂�A���̃V�~�����[�V�����𐧌䂷��B�ݒ�t�@�C���̋L�q�́A����̃p�[�g�ɕ�����Ă��邪�A�����ł�control���ɋL�q���ꂽ�p�����[�^�[���������Bcontrol���ȊO�Œ�`�����p�����[�^�[�ɂ��ẮAModule.f90���̋L�q�����Q�Ƃ��邱�ƁB�Ȃ��AFlag_�Ŏn�܂�ϐ��͘_���^�p�����[�^�[�Ȃ̂ŁA�u.true.�v���́u.false.�v�̂����ꂩ�����͂����B

Simulation_year       �V�~�����[�V����������N��
Flag_spinup_read  ���X�^�[�g���邩�ۂ��̃X�C�b�`
Flag_spinup_write ���X�^�[�g�t�@�C�����쐬���邩�ۂ��̃X�C�b�`
Flag_output_write �o�̓t�@�C�����쐬���邩�ۂ��̃X�C�b�`
Max_loc           ���z�ѕ��̈�Ђ̒���(m)
Depth             �y��w�ꖇ�̌���(mm)�B�y��͑S����30�w���݂���
STEP              �����̐����v�Z��A������DISK�ɕ�������ۂ̍ŏ��P�ʁB�ʏ��0.1�Ǝw��̂���(0.1m)�B
C_in_drymass      �o�C�I�}�X�̒P�ʊ����d�ʓ�����ɒY�f����߂銄���B
File_no           �o�̓t�@�C���̑��u�ԍ��B�����ԍ��𕡐���w�肵�Ȃ����ƁB�܂��A������Fortran�R���p�C���ł́A5�Ԃ��W�����́i�L�[�{�[�h�j�A6�Ԃ��W���o�́i��ʁj�Ɍ�������Ă���̂ŁA5��6�͎w�肵�Ȃ����ƁB
Fn_climate        �C�ۃf�[�^�̏��݁i�t���p�X�j
Fn_location       �S���̓y�n���Ɠy��f�[�^���t���p�X�Ŏw��B�v7�̃t�@�C�����t���p�X�Ŏw�肷��F���ʊC�m�}�X�N�E�W���E�y��A���x�h�E�ޏ�e���ʁE�ő�ې��ʁE��z���́E�ނ�_�B�e�f�[�^�t�@�C���̌`���́A��q�B
Fn_spnin          ���X�^�[�g���̓f�[�^���t���p�X�Ŏw��
Fn_spnout         ���X�^�[�g�o�̓f�[�^���t���p�X�Ŏw��


���������@SEIB-Viewer�@��������
SEIB-Viewer�́ASIEB-DGVM�̏o�͂��������邽�߂�Windows�p�A�v���P�[�V�����ł��B�C���X�g�[���͕K�v����܂��񂵁A���W�X�g�����g���܂���̂ŁA�t�H���_"SEIB-Viewer"��C�ӂ̏ꏊ�ɒu���āA���̒���SEIB-Viewer.exe�����s���邾���ŗ��p�ł��܂��B���̃A�v���P�[�V�����Ɋւ���A����ȏ�̐����́A�A�v���P�[�V��������Q�Ƃł���w���v�t�@�C���ɋL�q����Ă��܂��B

�Ȃ��ASEIB-Viewer��Visual Basic 2005 �ɂč쐬����Ă���܂��̂ŁA���̎��s�ɂ̓p�\�R���Ɂu.NET Framework 2.0�v���C���X�g�[������Ă���K�v������܂��B���̃C���X�g�[����Windows Update�𗘗p���Ă��\�ł����A���̃T�C�g�ihttp://msdn.microsoft.com/netframework/downloads/updates/default.aspx�j����.NET Framework Version 2.0 Redistributable Package����肷��̂���ԊȒP���Ǝv���܂��B�Ȃ��A32�r�b�g��:x86�A64�r�b�g��:x64�A64�r�b�g��:IA64�A�̂R��ނ̃p�b�P�[�W������܂��̂ŁA���g���̊��ɍ��킹���p�b�P�[�W����肵�Ă��������i�悭������Ȃ��ꍇ�ɂ́AWindows Update�𗘗p�����̂����S�ł��j�B
��Vista�A�܂��͂���ȏ�̃o�[�W������Windows�̏ꍇ�ɂ́A���̃C���X�g�[���͕K�v����܂���B


���������@main.f90���̃R�[�h�̗���@��������
main.f90�̓V�~�����[�^�[�̊�{���i�ł���Amain.f90�̈�ԊO���̃��[�v���P�V�~�����[�V�������ɑΉ����Ă���B���̂P���̃v���Z�X�ɂ����Ċe��f�ߒ��̃T�u���[�`���������Ăяo����邪�A�R�[�h�̉ǐ������߂邽�߂ɁASEIB-DGVM�ł́Amain.f90���Ăяo���ꂽ�T�u���[�`������A���̃T�u���[�`�����Ăяo�����Ƃ͋ɗ͔����Ă���B
�Ȃ����݂̃R�[�h�ł́A�w�肳�ꂽ�C�ۃf�[�^���J��Ԃ����͂����d�l�ƂȂ��Ă���B�Ⴆ��15�N���̋C�ۃf�[�^���w�肵���ꍇ�A15�N���ƂɋC�ۃf�[�^�������߂���ă��f���ɓ��͂����B

�������l�̐ݒ�i�v���O�����J�n���񂾂��Ăяo�����j
�ECall init_value
�e��ϐ��ɏ����l��^����

�ECall spinup_in
�K�v�ł���΃��X�^�[�g�t�@�C����ǂݍ���

��DAILY UPDATE OF FIELD PROPERTIES
�e��ϗʂ̃��Z�b�g�A�z��ϐ��ւ̓ǂݍ��݁A�ړ����ς̎Z�o

�����������W���[���Q�iDAILY�j
�ECall climate_stat
�C�ۊ��Ɋւ��鏔����z��ϐ��Ɏ�荞�݁A�܂�������p���Ĉړ����ς��Z�o����B

�ECall air
�C����O���ȂǁA��C�����Ɋւ��ϗʂ��Z�o

�ECall radiation
������PAR�ȂǁA�������Ɋւ��ϗʂ��Z�o

�ECall diffused_radiation
���z�ѕ��̉��������ɂ�����A�g�U���̑��΋��x���z�i���z�ѕ��̏�󂪊�j���Z�o����B

�ECall direct_radiation
�e�̂̊e�������C���[�ɂ����钼�ڌ��̑��΋��x�i���z�ѕ��̏�󂪊�j���Z�o����B�{�v���O�����ɂ����āA��������Ԍv�Z���Ԃ�����邽�߁A�����Ɉ��̂ݎ��s�����悤�ɂ��Ă���B�v�Z�͂ɗ]�͂�����̂Ȃ�΁A�������s���������ǂ����낤�B

�ECall floor_radiation
���ΏƓx�i�X�т̏�[�Ɣ�r�����j���A�ؖ{�̒蒅Cell���ƂɎZ�o����B

�ECall ir_index
�ؖ{�w�S�̂��A�܂����{�w�S�̂��A���z������������������̂��������W���iir_tree�Air_grass�j���Z�o����B�܂��A���{�w�̍ŏ㕔�ɂ�����������L�����˗ʁipar_grass�j���A�����ŎZ�o�����

�ECall fire_regime
�X�щ΍Ѓ��W���[���B�A�t���J�嗤�̏ꍇ�ɂ�fire_regime2��Call�����B

�������ߒ��iDAILY�j
�ECall photosynthesis_condition
���������x�̎Z�o�ɕK�v�ȏ��ϐ��i���O�a���̌��������x�Ȃǁj���Z�o����

�ECall photosynthesis
1���̌������ʂ��Z�o����i���{�E�ؖ{�Ƃ��j

�ECall lai_optimum
���{���C���[�ɂ�����A�œK�t�ʐώw���̎Z�o

�ECall maintenance_resp
�ێ��ċz�i���{�E�ؖ{�Ƃ��j

�ECall turnover
�^�[���I�[�o�[�i���{�E�ؖ{�Ƃ��j

�ECall leaf_season
�t�̃t�F�m���W�[�̐؂�ւ��ƁA���̐؂�ւ��̑O��ɂ�����A�e�픽���i���{�E�ؖ{�Ƃ��j

�ECall growth_wood
�ؖ{�̐����̂����A1���x�[�X�ŃV�~�����[�g���镔���i�W�t�A�׍��̐����A�ɐB�j

�ECall growth_grass
���{�w�̐����E�ɐB

�ECall spacial_limitation
�e�ؖ{�̂̎��͂̋�ԓI�]�T���Z�o����v���O�����B�ȉ��̃T�u���[�`��growth_trunc�ɂ����āA�����̐�������������g�傳����O�ɁA�Ăяo���B

�ECall growth_trunc
���̔�听���ƁA����ɔ����������Ǝ������̐���

�ECall decomposition
�y��L�@���̕������W���[��

�������ߒ��iANNUALY�j
�ECall mortality
�X�̖ؖ{�ɂ��āA���S�𔻒f����

�ECall pft_present
���݁A���̉��z�ѕ��ɑ��݂���PFT�����X�g�A�b�v����

�ECall crown_adjust
�X�̖ؖ{�ɂ��āA�҂��̈����t�Q���C���[���͂�グ������

�ECall crown_shake
�X�̖ؖ{�ɂ��āA���̎����̐��������̈ʒu���A����Ԃ̋󂢂Ă�������Ɉړ�������

�ECall ground_vacant
���͂̋�ԓI�]�T�Ȃǂ����ɁA�蒅�\�ȃO�����h���b�V�������X�g�A�b�v����B���L�̃T�u���[�`��establish�ŁA�t���̒蒅���s�킹��O�ɌĂяo�����

�ECall establish
�t���̒蒅

�ECall crown_coverage
�ؖ{�̎����ɂ��n�\�ʂ̃J�o�[�����Z�o

�ECall grass_priority
���N�ɗD��^�Ƃ��鑐�{PFT���A�O�N�x�̒P�ʖʐϓ�����NPP���r���邱�ƂŌ��߂�

�ECall biome_determine
LAI��D��PFT�Ƃ������������ɁA���݂̐��Ԍn�^�C�v�����肷��

��DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS
�ECall net_radiation
���ˎ��x��A���x�h�̌v�Z

�ECall waterbudget
�����x�̌v�Z

��Daily increment of litter and fuel pools
��{�I�ɁA���^�[�v�[����Fuel�v�[���̑����́A���C���ȉ��̊e�T�u���[�`���ŎZ�o���ꂽdaily litter fluxes�����ɁA���̃p�[�g�ŏW�v�����B���̃��f���̓��^�[�v�[����Fuel�v�[�����R�����A���̍\�������X��₱�������߁A�܂Ƃ߂ďW�v���ăR�[�h�̍���������邱�Ƃ��A���̏W�v���@�̖ړI�B

��MAKE OUTPUT FILES
�ECall output_*
�V�~�����[�V�������ʂ̏o��


���������@��̓v���O�����@��������
���X�ѓ��Ԃ̊�{��� (analysis.f90)
analysis.f90�́A�o�̓t�@�C���̈�ł���output_forest.txt��p���āA�X�ѓ��ԂɊւ����{�����Z�o����B�Z�o�������́A�T�C�Y�N���X���̖ؖ{���x�A�T�C�Y�N���X���̊���呬�x�A�T�C�Y�N���X���̎��S���̂R�ł���B�Ȃ��A�����ŃT�C�Y�N���X�Ƃ́A���̒��a�����Ƃ����ؖ{�T�C�Y�̎w�W�ł���B
�v���O��������Max_loc�AMax_no�ASimulation_year�ɂ̓V�~�����[�V�����̐ݒ�t�@�C���Ɠ����l����͂���BInterval1�ɂ́A���S���Ɛ������x�̎Z�o�ɍۂ��ĐX�т��r����Ԋu��'�N'�œ��͂���B�Œ�l��1�ł��邪�A�����Ȗʐςł̃V�~�����[�V�����̏ꍇ�A���̊Ԋu���]��ɒZ���ƌ덷���傫���Ȃ�BInterval2�ɂ́A��͂̌J��Ԃ��Ԋu�ł���i�N�j�B�Ⴆ�΁A�����50�Ɛݒ肵�A1000�N���̃V�~�����[�V�����f�[�^�ɂ��Ď��s�����ꍇ�A�V�~�����[�V�����J�n��50, 100, 150, ...., 1000�N�̌v20��ɂ��ĉ�͂��s���A���̕��ϒl���o�͂���BSizeclass_num�ł̓T�C�Y�N���X�̌����`���ASize_set�ł́A�T�C�Y�N���X�́u�؂���v���`����B�Ⴆ�΁A���a��10-15, 15-20, 20-25, 25�ȏ��4�N���X��ݒ肷��ꍇ�ɂ�Sizeclass_num�́u4�v�ASize_set�́u10,15,20,25�v�ƒ�`����B


�����z�X�т̉����iViewer_snapshot.pov, Viewer_animation.pov�j
�o�̓t�@�C���̈�ł���output_forest.txt��p���ĉ��z�ѕ��̎O�����\����`�悷�邱�Ƃ��ł���B�܂�pov-ray�ihttp://www.povray.org/�j�Ƃ��������̃����_�����O�\�t�g���C���X�g�[������B����output_forest.txt�Ɠ����f�B���N�g���ɓY�t��Viewer_snapshot.pov�������A���s������B������#declare year���ŕ\�����������N�i�V�~�����[�V�����J�n��̔N���j����͂��A����̃v���_�E�����j���[�ŏo�͉摜�̑傫�������w�肵�A�uRun�v�Ə����ꂽ�{�^���������Ή��z�ѕ������������B���̉摜�́A�����f�B���N�g���Ƀr�b�g�}�b�v�`���Ŏ����ۑ�����Ă���B

Viewer_animation.pov�͓���쐬�p�̃R�[�h�ŁA�P�N���Ƃ̘A�������X�ё��ς��o�͂���B�Ⴆ��100�N�Ԃ̐X�ѕϓ���`������̂ł���΁A�v���O����14�s�ڂ��u#declare Flame = 100;�v�Ə����A���s�I�v�V�����̓��͗��i��ʉE��̋󗓁j�Ɂu+kff100�v�Ɠ��͂��A�uRun�v�{�^���������B����ƘA�������r�b�g�}�b�v�t�@�C�������X�Əo�͂����̂ŁA�����AviMaker�ihttp://www.vector.co.jp/soft/win95/art/se121264.html�j�Ȃǂ�p���ē���֕ϊ�����B

pov-ray�Ɋւ��ē��{��ŏ����ꂽ���́A�ȉ��̃T�C�g�Ȃǂɂ���B
http://www.geocities.co.jp/SiliconValley/5518/
http://nishimulabo.edhs.ynu.ac.jp/~povray/3.5jp/index.htm
http://www.arch.oita-u.ac.jp/a-kse/povjp/index.htm


���������@�G�f�B�^�ɂ��ā@��������
Windows��ł̃R�[�h�̕\����ҏW�ɂ́u�G�ۃG�f�B�^�v�Ƃ����e�L�X�g�G�f�B�^���֗��ł��B�ݒ胁�j���[����܂�Ԃ���������100�ȏ�ɐݒ肵�AFortran�p������`��ǂݍ��܂���ƁA�Y�킩������₷���\������܂��B������`�t�@�C����ǂݍ��܂�����@�́A�u���̑����t�@�C���^�C�v�ʂ̐ݒ聨�����\�����ǂݍ��݁v�ƃ��j���[�I��������A���̋����\���t�@�C�����w�肵�܂��B

�u�G�ۃG�f�B�^�v�̓V�F�A�E�F�A�ł��̂Ōp�����ė��p����ۂɂ͑������K�v�ł����A�w���̕��́A�w�Ђɂ���Ԃ̂ݑ������Ə������u�A�J�f�~�b�N�t���[���x�v�����݂��܂��B���̐��x�𗘗p�����ۂɂ́A�ȉ��̃T�C�g�́u�T�|�[�g�v�y�[�W�ŕK�v��������͂��Đ\������ł��������B

�G�ۃG�f�B�^�T�C�g�i�R�[�h�͂����炩��j
http://hide.maruo.co.jp/

�����\����`�͂����炩��
http://hide.maruo.co.jp/lib/hilight/index.html


���������@���s�R�[�h�ɂ��ā@��������
readme_j.txt�i���̕��́j�Areadme_e.txt�A������POV-RAY�v���O�����������āA�S�Ẵt�@�C���́A���s�R�[�h��UNIX�W���ł���uLF�`���v���g�p���Ă��܂��BWindows��Ńt�@�C�����K�؂ɕ\���ł��Ȃ��ۂɂ́A�܂��͉��s�R�[�h��Windows�W���ł���uCR+LF�`���v�ɕϊ����Ă݂Ă��������B���s�R�[�h�̕ϊ��ɂ́A���̃t���[�E�F�A���֗��ł��ihttp://www.vector.co.jp/soft/win95/util/se296459.html�j�B���l�ɁAUNIX�T�[�o�[��ȂǂŁA�t�@�C�����������\������Ȃ��ۂ��A�܂��͉��s�R�[�h�����m���߉������B


���������@���̑��̒��Ӂ@��������
 ���z�z�����R�[�h�͈�n�_�v�Z�p�ł����Astart_point.f90�����������邱�Ƃɂ�蕡���n�_�̘A���V�~�����[�V�������\�ɂȂ�܂��B�A���A����Ȃ�ɋ��͂Ȍv�Z�����K�v�ƂȂ�܂��B�����n�_�v�Z�p��start_point.f90�A�܂��X�p�R���փW���u�𓊓�����ۂ̊e��t�@�C��(makefile����s�o�b�`)���K�v�ȏꍇ�͂��A�������� 
