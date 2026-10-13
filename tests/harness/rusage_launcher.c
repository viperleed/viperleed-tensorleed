/* Run a program and report its resource usage.
 *
 * Usage: rusage_launcher RESULT_FILE PROGRAM [ARGS...]
 *
 * Forks, execs PROGRAM (inheriting stdin/stdout/stderr and cwd), waits for
 * it with wait4() and writes "<wait status> <user s> <sys s> <maxrss kB>"
 * to RESULT_FILE. Needed because on Linux ru_maxrss survives execve: a
 * program started directly from Python would report at least the Python
 * process's peak RSS. This launcher is small, so its own footprint (~1 MB)
 * is the floor of the measurement. The child is killed if the launcher dies
 * (timeout handling in harness/runner.py).
 */
#include <signal.h>
#include <stdio.h>
#include <sys/prctl.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s RESULT_FILE PROGRAM [ARGS...]\n", argv[0]);
        return 125;
    }
    pid_t pid = fork();
    if (pid < 0) {
        perror("rusage_launcher: fork");
        return 125;
    }
    if (pid == 0) {
        prctl(PR_SET_PDEATHSIG, SIGKILL);
        execv(argv[2], argv + 2);
        perror("rusage_launcher: execv");
        _exit(127);
    }
    int status;
    struct rusage ru;
    if (wait4(pid, &status, 0, &ru) < 0) {
        perror("rusage_launcher: wait4");
        return 125;
    }
    FILE *f = fopen(argv[1], "w");
    if (!f) {
        perror("rusage_launcher: result file");
        return 125;
    }
    fprintf(f, "%d %ld.%06ld %ld.%06ld %ld\n", status,
            (long)ru.ru_utime.tv_sec, (long)ru.ru_utime.tv_usec,
            (long)ru.ru_stime.tv_sec, (long)ru.ru_stime.tv_usec, ru.ru_maxrss);
    fclose(f);
    return 0;
}
