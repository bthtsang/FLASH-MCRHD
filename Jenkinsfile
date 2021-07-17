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
    
    stage('CartRadDiff'){
        dir('FLASH4'){
            sh 'python3 bin/setup.py CartRadDiff -auto -3d -maxblocks=200 -debug -objdir=cartraddiff'
            dir('cartraddiff'){
	        stage('Build'){
                    sh 'cp /home/jenkins/Makefile.h .'
                    sh 'make -j'
		}
		stage('Test'){
                    sh 'mpirun -np 8 ./flash4'
                    sh 'python3 urad_from_mcps.py cartimc_1M_hdf5_chk_ 1 3'
                    archiveArtifacts artifacts: 'FLASH4/cartraddiff/mcp_urad.pdf'
		}
            }
        }
    }
}
