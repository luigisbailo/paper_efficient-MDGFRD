// author luigisbailo


void printPos_per_nl ( struct particle_nl *particles, int *partList, int N){

    printf("Part.\tTime\tTau\tShell\tx\ty\tz\tx_ex\ty_ex\tz_ex\tper_x\tper_y\tper_z\tRadius\tBURST\n");

    for (int count=0; count<N; count++){

        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf\t%d\t%d\n",
               particles[partList[count]].label,particles[partList[count]].time,
               particles[partList[count]].tau_exit,particles[partList[count]].shell,
               particles[partList[count]].pos[0], particles[partList[count]].pos[1],
               particles[partList[count]].pos[2], particles[partList[count]].pos_exit[0],
               particles[partList[count]].pos_exit[1], particles[partList[count]].pos_exit[2],
               particles[partList[count]].pos_period[0],particles[partList[count]].pos_period[1],
               particles[partList[count]].pos_period[2],
               particles[partList[count]].radius, particles[partList[count]].burst,
               particles[partList[count]].gf);

    }

}

