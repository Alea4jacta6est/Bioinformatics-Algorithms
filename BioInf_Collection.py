def Needleman_Wunsch(seq_1, seq_2, indel=-1, missmatch=-1, match=1, cost_mat=None):
    import numpy as np
    import itertools

    def cost_func(el_1, el_2, indel=indel, missmatch=missmatch, match=match, cost_mat=cost_mat):
        if cost_mat:
            return cost_mat.acces(el_1, el_2)
        elif el_1=='-' or el_2=='-':
            return indel
        else:
            return match if el_1==el_2 else missmatch
    
    Dist_mat=np.empty([len(seq_1)+1, len(seq_2)+1])
    Path_mat=np.empty([len(seq_1)+1, len(seq_2)+1])
    
    Dist_mat[:,0]=[x*indel for x in range(0, len(seq_1)+1)]
    Dist_mat[0,:]=[x*indel for x in range(0, len(seq_2)+1)]
    
    Path_mat[:,0]=[0]*(len(seq_1)+1)
    Path_mat[0,:]=[1]*(len(seq_2)+1)

    for i, j in itertools.product(range(1, len(seq_1)+1), range(1, len(seq_2)+1)):
    
        Dist_mat[i, j]=max(Dist_mat[i-1,j]+cost_func(seq_1[i-1], '-'), Dist_mat[i, j-1]+cost_func('-', seq_2[j-1]), 
                           Dist_mat[i-1,j-1]+cost_func(seq_1[i-1], seq_2[j-1]))
        
        Path_mat[i, j]=np.argmax([Dist_mat[i-1,j]+cost_func(seq_1[i-1], '-'), Dist_mat[i, j-1]+cost_func('-', seq_2[j-1]), 
                           Dist_mat[i-1,j-1]+cost_func(seq_1[i-1], seq_2[j-1])])

    Alignment=[[],[]]
    
    i,j=len(seq_1),len(seq_2)
    
    while (i,j)!=(0,0):
        
        if Path_mat[i, j]==0: #Deletion
            Alignment[0].append(seq_1[i-1])
            Alignment[1].append('-')
            i-=1
        elif Path_mat[i, j]==1: #Insertion
            Alignment[0].append('-')
            Alignment[1].append(seq_2[j-1])
            j-=1
        else: #Match/Missmatch
            Alignment[0].append(seq_1[i-1])
            Alignment[1].append(seq_2[j-1])
            i-=1
            j-=1
            
        
    Alignment[0]=Alignment[0][::-1]
    Alignment[1]=Alignment[1][::-1]
    return Dist_mat, Alignment