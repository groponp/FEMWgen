#!-----------------------------------------------------------------------------------------------#!
#! FEMWgen  : A computational tool to generate multiple windows                                  #!
#!            for free-energy calculations.                                                      #!
#!                                                                                               #!
#! @uthors : Alexandre Suman de Araujo, Luiz Fernando Zonetti e Ropón-Palacios G.                #!
#! E-mail  : <asaraujo@ibilce.unesp.br>.                                                         #!
#! Versão  : 1.0b                                                                                #!
#! Data    : Tue 8 Nov, 2022.                                                                    #!
#!-----------------------------------------------------------------------------------------------#!

#! Observações:
#!--------------------------------
#! 1) Esse código foi criado a partir de um conjunto de scripts que o Alexandre montou em 2010.
#! 2) Agora iremos tentar criar um único script que executa todas as funções de forma geral dada uma membrana e um peptideo.

#! ChangeLogs:
#!--------------------------------
#! 1) 0.1a - Nesse script alteramos o scritp roda_abf para somente criar as janelas do ABF. Além do mecanismo de empurrar a proteina pra dentro da membrana
#! 2) implementamos um metodo de "abrir" a membrana para acomodar a proteina no interior da mesma de forma correta.
#! 3) 0.3a - Esse script consistia de um conjunto de arquivos e não era geral para qualquer membrana e peptideo.
#! 4) 0.4a - Nesse iremos tentar generalizar e compactar o 0.3a.
#! 5) 0.5a - Nesse foi adicionada command line, para ser mais facil usar con VMD. 22:36 pm, Tue 8, 20.22.
#! 6) 0.5a - Foi trocado `echo` pelo `puts` 23:51 pm, Tue 8, 2022. 
#! 8) 1.0b - Foi adicionado a completo command line procedure, para poder usar en BATCH mode. 14:35 pm, Wed 9 Nov, 2022. 

#! Notas:
#!--------------------------------
#! 27/10/2014: 
#! 1) O arquivo em_janelas.conf é necessário no processo de encolhimento da membrana e deve ser copiado e adaptado ao seu sistema junto com os outros arquivos de simulação.
#! 2) Aqui estamos considerando que nosso sistema comeca com o peptideo na parte NEGATIVA e vai para a positiva.
#! 3) Outro detalhe é que o .pdb deve conter primeiro a proteína e a membrana e depois a água e íons para não dar problema na hora de criar as restrições.


package provide FEMWgen 1.0 

#! Add command line procedure
#!-----------------------------
#! set env(VMDNOMANGLEATOMNAMES) 1

proc femwgen_usage { } {
   puts "#!-------------------------------------------------------------------------------------------------#!"
   puts "#! FEMWgen v1.0b : A Computational Tool to Generate Multiple Windows for Free-energy Calculations. #!"
   puts "#! \[INFO   \] usage: source FEMWgen.tcl.                                                            #!"
   puts "#! \[INFO   \] usage: FEMWgen -psf<PSF> -pdb<PDB> -n<N> -coj<CPJ> -bin<BIN> -lresn<LRESN>...         #!"      
   puts "#! \[INFO   \] @uthors : Suman de Araujo A., Fernando Zonnetti and Ropón-Palacios G.                 #!"
   puts "#! \[INFO   \] contact to: <asaraujo@ibilce.unesp.br>.                                               #!"
   puts "#!-------------------------------------------------------------------------------------------------#!"
   puts "Options:"
   puts "    -psf<PSF>         File containing topology. Default 'sistema.psf'."
   puts "    -pdb<PDB>         File containing coordinates. Default 'sistema.pdb'."
   puts "    -n<N>.            Number of windows. Default 60."
   puts "    -cpj<CPJ>         Center of first windows. Default -38."
   puts "    -bin<BIN>         Spacing for each windows. Default 2."
   puts "    -lresn<LRESN>     Resname lipid component into membrane. Default 'POPC POPG'."
   puts "    -pref<PREF>       Reference atoms from protein. Default 'resname TRP and name CH2'."
   puts "    -lref<LREF>       Reference atoms from lipids. Default 'resname POPC POPG and name P'."
   puts "    -lrst<LRST>       Serial key of atoms to restraint during equilibration. Default 'serial 38224 40079 39148 39811 20055 14892 19787 18877'."
   puts "    -charm++<CHARM++> Boolean option to use o not charmrun mode with running NAMD. Default 'false'."
   puts "    -nprocs<nprocs>   Number of CPUs that will be use into energy minimization. Default 2."
   puts "    -h                Print it message."
   error ""

}

proc FEMWgen { args } {

   set cmdlinelength [llength $args]


   if { $cmdlinelength ==! 11 || [lindex $args 0] == "-h" } {

      femwgen_usage
   }

   #! Set the defaults.
   #!-----------------------
   set  psf_file        "sistema.psf"                                               ;#! Arquivo de topología.
   set  pdb_file        "sistema.pdb"                                               ;#! Arquivo de coordenadas.
   set  N                60                                                         ;#! Número de janelas.
   set  CPJ              -38                                                        ;#! Centro da primera janela.
   set  BIN              2                                                          ;#! Largura da janela.
   set  lipids_in_memb   "POPC POPG"                                                ;#! Resname dos lípidos que componem na membrana.
   set  sel_ref_prot     "resname TRP and name CH2"                                 ;#! Seleção de referencia na proteina.
   set  sel_ref_memb     "resname POPC POPG and name P"                             ;#! Seleção de referencia na membrana.
   set  sel_rest_memb    "serial 38224 40079 39148 39811 20055 14892 19787 18877"   ;#! Seleção para restringir o movimiento da membrana no porceso de equi.
   set  charmrun         "false"                                                    ;#! Neste opcao so é cuando vc tem namd2 compilado con charmrun.
   set  nprocs           2                                                          ;#! É número de CPUs pra usar na minimizacao.


   #! Parse options.
   #!--------------------
   for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {

      set arg [ lindex $args $argnum ]              #! Get command.
      set val [ lindex $args [expr $argnum + 1] ]   #! Get value to this command.
      switch -- $arg {

         "-psf"        { set psf_file $val; incr argnum; continue }
         "-pdb"        { set pdb_file $val; incr argnum; continue }
         "-n"          { set N $val; incr argnum; continue }
         "-cpj"        { set CPJ $val; incr argnum;continue }
         "-bin"        { set BIN $val; incr argnum; continue }
         "-lresn"      { set lipids_in_memb $val; incr argnum; continue }
         "-pref"       { set sel_ref_prot $val; incr argnum; continue }
         "-lref"       { set sel_ref_memb $val; incr argnum; continue }
         "-lrst"       { set sel_rest_memb $val, incr argnum; continue }
         "-charm++"    { set charmrun $val; incr argnum; continue }
         "-nprocs"     { set nprocs $val; incr argnum }

         default { puts "unknow option: $arg"; return }

      }
   }

proc writeFiles { sel ofile  mode} {

   if { $mode == "both" } {

      package require psfgen
      resetpsf 
      $sel writepsf $ofile
      $sel writepdb $ofile 

   } elseif { $mode == "psf" } {

      package require psfgen
      resetpsf
      $sel writepsf $ofile 

   } elseif { $mode == "pdb" } {

      $sel writepsf $ofile 

   } else {

      puts "\[INFO   \] Not mode select."
   }


}

#! Inicia o loop de geração de janelas.
#!----------------------------------------
for {set i 1} {$i<=$N} {incr i} {

   #! Carrega o sistema original e posiciona o peptideo no centro da primeira janela.
   #!---------------------------------------------------------------------------------
   mol delete top
   mol new $psf_file
   mol addfile $pdb_file waitfor all
   set prot [atomselect top protein]
   set ref_memb [atomselect top $sel_ref_memb]
   set ref_prot [atomselect top $sel_ref_prot]
   set dz [expr ($CPJ - [lindex [measure center $ref_memb] 2]) - ([lindex [measure center $ref_prot] 2] - [lindex [measure center $ref_memb] 2])]
   lappend moveprot 0 0 $dz
   $prot moveby $moveprot
   unset moveprot   
   
   #! Gerando as pastas e copiando os arquivos necessários e entra na pasta.
   #!--------------------------------------------------------------------------
   mkdir janela_$i

   #! PONTO DE ALTERACAO.
   #!--------------------------
   set files [glob configs/*]
   foreach file $files{

      cp $file janela_$i/
   }
   #cp em_janela.conf em.conf md_eq.conf md_abf.conf md_abf.in par_all27_prot_na.prm par_all36_lipid.prm janela_$i/
   #cp configs/* janela_$i/
   cd janela_$i

   #! Translate o peptídeo.
   #!--------------------------
   set CJ [expr $CPJ+($i-1)*$BIN]
   set dz [expr $CJ-$CPJ]
   set prot [atomselect top protein]
   lappend moveprot 0 0 $dz
   $prot moveby $moveprot
   unset moveprot

   #! Cria a seleção sem águas a 2 A da proteina e salva os arquivos .pdb .psf.
   #!----------------------------------------------------------------------------
   set saida [atomselect top "not same residue as water within 2 of protein"]
   [writeFiles $saida $pdb_file "both"]

   #! Remove o sistema atualmente carregado.
   #!---------------------------------------------
   mol delete top
   
   #! Carrega novo sistema gerado e cria selecoes novamente.
   #!------------------------------------------------------------
   mol new $psf_file
   mol addfile $pdb_file waitfor all
   set prot [atomselect top protein]

   #! Verifica se está dentro da membrana.
   #! Pega z minimo e maximo da proteina.
   #!---------------------------------------------
   set prot_min [lindex [lindex [measure minmax $prot] 0] 2]
   set prot_max [lindex [lindex [measure minmax $prot] 1] 2]

   #! Seleciona os extremos de cada lado da membrana e define sua posicao media em Z.
   #!-----------------------------------------------------------------------------------
   set p_menos [atomselect top "resname $lipids_in_memb and name P and z < 0"]
   set p_mais [atomselect top "resname $lipids_in_memb and name P and z > 0"]
   set p_menos_z [lindex [measure center $p_menos] 2]
   set p_mais_z [lindex [measure center $p_mais] 2]

   #! Calcula as distancias entre a ponta e cauda da prot e os extremos da membrana.
   #!-------------------------------------------------------------------------------------
   set z1 [expr $p_menos_z-$prot_max]
   set z2 [expr $p_mais_z-$prot_min]

   #! Se a proteina estiver na regiao interior da membrana realiza o processo de abertura, se não, cria as restricoes, 
   #! atualiza .in e vai para proxima janela.
   #!-------------------------------------------------------------------------------------------------------------------
   if {$z1 <= 2 && $z2 >= -2} {

      puts "\[INFO   \] Entrou na abertura da membrana."

      #! Cria a selecao dos atomos de referencia e pega seu valor X.
      #!----------------------------------------------------------------
      set atomos_ref [atomselect top $sel_rest_memb]
      set pos_x_ref [$atomos_ref get x]

      #! Pega o centro e a largura em X da proteina.
      #!-------------------------------------------------
      set centro_x_prot [lindex [measure center $prot] 0]
      set largura_x_prot [expr [lindex [lindex [measure minmax $prot] 1] 0] - [lindex [lindex [measure minmax $prot] 0] 0]]

      
      #! Define as bandas para a abertura da membrana.
      #!-------------------------------------------------
      puts "\[INFO   \] Definindo as bandas."
      
      #! Seleciona a membrana.
      #!--------------------------------
      set memb [atomselect top "resname $lipids_in_memb"]

      #! Corre cada elemento da selecao da membrana e verifica a que banda pertence.
      #!-------------------------------------------------------------------------------
      set residue_mol [lindex [$memb get residue] 0]
      if {[lindex [measure center [atomselect top "residue $residue_mol"]] 0] < $centro_x_prot} {
         lappend residues_banda_inf $residue_mol
      } else {
         lappend residues_banda_sup $residue_mol
      }

      foreach temp [$memb get residue] {
         if {$temp != $residue_mol} {
            set residue_mol $temp
            if {[lindex [measure center [atomselect top "residue $residue_mol"]] 0] < $centro_x_prot} {
               lappend residues_banda_inf $residue_mol
            } else {
               lappend residues_banda_sup $residue_mol
            }
         }
      }
      
      set banda_sup [atomselect top "residue $residues_banda_sup"]
      set banda_inf [atomselect top "residue $residues_banda_inf"]

      #! Movimenta as bandas.
      #!-----------------------------
      puts "\[INFO   \] Abrindo a membrana."

      lappend d_sup [expr $largura_x_prot/2] 0 0
      lappend d_inf [expr (-1)*$largura_x_prot/2] 0 0
      $banda_sup moveby $d_sup
      $banda_inf moveby $d_inf
      unset d_sup
      unset d_inf
      
      #! Seleciona tudo e escreve o pdb e o freeze.
      #!-----------------------------------------------
      set tudo [atomselect top all]
      [writeFiles $tudo $pdb_file "pdb"]
      $tudo set beta 0
      set freezing [atomselect top protein]
      $freezing set beta 1
      [writeFiles $tudo "myfixedatoms.pdb" "pdb"]

      #! Processo de fechamento da membrana.
      #!--------------------------------------------
      puts "\[INFO   \] Comeca a fechar a membrana."
      
      #! Cria a pasta temporaria e copia tudo para dentro.
      #!-----------------------------------------------------
      mkdir temp
      #! PONTO DE ALTERACAO.
      #!----------------------------
      cp $psf_file $pdb_file myfixedatoms.pdb em_janela.conf temp/
      cd temp

      #! Inicia o loop de encolhimento.
      #-----------------------------------
      set gatilho 1
      set rmsd_ant 10000.0
      while {$gatilho} {
         puts "\[INFO   \] Entrou no loop de fechamento."
            
         #! Minimiza.
         #!--------------------------------
         if { $charmrun == "true" } {

            charmrun +p$nprocs namd2 em_janela.conf >& min.log 
            puts "\[INFO   \] Olhe min.log file to check out progresso."

         } else {

            namd2 +p$nprocs em_janela.conf >& min.log 
            puts "\[INFO   \] Olhe min.log file to check out progresso."
         }


         #! Remove a molecula carregada e carrega a saida da minimizacao.
         #!---------------------------------------------------------------
         mol delete top
         mol new $psf_file
         #! ATENÇÃO.
         #!---------------------
         mol addfile em.coor waitfor all

         #! Recria a selecao de tudo.
         #!-------------------------------
         set tudo [atomselect top all]
         
         #! Recria a selecao dos atomos de referencia e pega seu valor X.
         #!----------------------------------------------------------------
         set atomos_ref [atomselect top $sel_rest_memb]
         set pos_x_temp [$atomos_ref get x]

         #! Calcula o RMSD entre as distancias de referencia e as atuais.
         #!----------------------------------------------------------------
         set rmsd 0
         for { set j 0 } { $j < [llength $pos_x_ref] } { incr j } {
            set rmsd [expr $rmsd + (([lindex $pos_x_temp $j]-[lindex $pos_x_ref $j])*([lindex $pos_x_temp $j]-[lindex $pos_x_ref $j]))]
         }
         set rmsd [expr sqrt($rmsd/[llength $pos_x_ref])]

         if {$rmsd > $rmsd_ant} {
            set gatilho 0
            cp sistema_ant.pdb $pdb_file
            exec /bin/sh -c "echo $rmsd >> rmsd.dat"
         } else {
            set rmsd_ant $rmsd
            exec /bin/sh -c "echo $rmsd >> rmsd.dat"
         }
            
            #! Salva o sistema minimizado.
            #!--------------------------------
            [writeFiles $tudo "sistema_ant.pdb" "pdb"]

            #! Cria novamente as bandas que serão movidas.
            #!--------------------------------------------
            set banda_sup [atomselect top "residue $residues_banda_sup"]
            set banda_inf [atomselect top "residue $residues_banda_inf"]

            #! Move as bandas.
            #!-----------------------
            $banda_sup moveby {-0.2 0 0}
            $banda_inf moveby {0.2 0 0}
            
            #! Escreve o .pdb da membrana encolhida e cria o arquivo de freeze.
            #!--------------------------------------------------------------------
            [writeFiles $tudo $pdb_file "pdb"]
            $tudo set beta 0
            set freezing [atomselect top protein]
            $freezing set beta 1
            [writeFiles $tudo "myfixedatoms.pdb" "pdb"]
            
         }
         
      
      #! Copia o sistema.pdb para a pasta anterior, desce uma pasta e apaga o temp.
      #!-----------------------------------------------------------------------------
      cp ../$pdb_file ../sistema_orig.pdb
      cp $pdb_file rmsd.dat ../
      cd ../
      rm -rf temp/
   }
      
   puts "\[INFO   \] Passou pela verificacao da abertura da membrana."
   #! Remove o sistema atualmente carregado.
   #!-------------------------------------------
   mol delete top
   
   #! Carrega a janela gerada.
   #!-------------------------------
   mol new $psf_file
   mol addfile $pdb_file waitfor all

   #! Cria arquivos de restraint e freezing.
   #! Cria arquivo myfixedatoms.pdb de freezing.
   #!--------------------------------------------
   set tudo [atomselect top all]
   $tudo set beta 0
   set freezing [atomselect top protein]
   $freezing set beta 1
   [writeFiles $tudo "myfixedatoms.pdb" "pdb"]

   #! Cria arquivo restraint.pdb de restraint dos fosforos da membrana.
   #!--------------------------------------------------------------------
   $tudo set beta 0
   set rest [atomselect top $sel_rest_memb]
   $rest set beta 1
   [writeFiles $tudo "restraint.pdb" "pdb"]      
   
   #! Atualiza o arquivo md_abf.in.
   #! Calcula o inicio e fim da janela.
   #!--------------------------------------
   set ij [expr $CJ-$BIN/2]
   set fj [expr $CJ+$BIN/2]

   #! Pega serial dos átomos de referência para o ABF.
   #!-----------------------------------------------------
   set main_serial [ [atomselect top $sel_ref_prot] get serial]
   set ref_serial [ [atomselect top $sel_ref_memb] get serial]

   #! Substitui tudo no arquivo .in.
   #!-----------------------------------
   exec /bin/sh -c "cat md_abf.in|sed s/SS/'$fj'/ > temp1.in"
   exec /bin/sh -c "cat temp1.in|sed s/II/'$ij'/ > temp2.in"
   exec /bin/sh -c "cat temp2.in|sed s/MMM/'$main_serial'/ > temp3.in"
   exec /bin/sh -c "cat temp3.in|sed s/RRR/'$ref_serial'/ > md_abf.in"
   exec /bin/sh -c "rm -rf temp*.in"

   #! Volta para a pasta original.
   #!---------------------------------
   cd ../

}
}

#Sai do VMD
#exit