node{
    
    stage('Organize'){
        sh "rm -rf ./*"
        
        dir('FLASH-MCRHD'){
	    checkout scm
            sh 'ls'
        }
        
        dir('FLASH4'){
            git branch: 'master', credentialsId: 'github_personal_access_token', url: 'https://github.com/srichers/FLASH4.6.2.git'
            sh 'ls'
            sh "cp -r ../FLASH-MCRHD source/physics/MonteCarloRHD"
            sh "cp -r ../FLASH-MCRHD/TestDirectories/* source/Simulation/SimulationMain"
            sh "ls source/physics/MonteCarloRHD"
            sh "ls source/physics/MonteCarloRHD/MonteCarloRHDMain"
        }
        
    }
    
    stage('RadiativeShock'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH4'){
            sh 'python3 bin/setup.py RadiativeShock -auto -1d -maxblocks=100 -debug -objdir=radiativeshock1d'
            dir('radiativeshock1d'){
                sh 'cp /home/jenkins/Makefile.h .'
                sh 'make -j'
                sh 'mpirun -np 8 ./flash4'
		sh 'python3 radshock_super_1d.py'
		archiveArtifacts artifacts: '*.pdf'
            }
        }
    }}

    stage('CartRadEqm'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH4'){
            sh 'python3 bin/setup.py CartRadEqm -auto -3d -maxblocks=200 -debug -objdir=cartradeqm'
            dir('cartradeqm'){
                sh 'cp /home/jenkins/Makefile.h .'
                sh 'make -j'
                sh 'mpirun -np 8 ./flash4'
                sh 'python3 plot_temperature.py'
                archiveArtifacts artifacts: 'cartradeqm_TfluidTrad.pdf'
            }
        }
    }}

    stage('CartRadDiff'){catchError(buildResult: 'UNSTABLE', stageResult: 'FAILURE'){
        dir('FLASH4'){
            sh 'python3 bin/setup.py CartRadDiff -auto -3d -maxblocks=200 -debug -objdir=cartraddiff'
            dir('cartraddiff'){
                sh 'cp /home/jenkins/Makefile.h .'
                sh 'make -j'
                sh 'mpirun -np 8 ./flash4'
                sh 'python3 urad_from_mcps.py cartimc_1M_hdf5_chk_ 1 3'
                archiveArtifacts artifacts: 'mcp_urad.pdf'
            }
        }
    }}
}
