" Vimball Archiver by Charles E. Campbell, Jr., Ph.D.
UseVimball
finish
autoload/pyclewn/start.vim	[[[1
223
" vi:set ts=8 sts=4 sw=4 et tw=80:
" Pyclewn run time file.
" Maintainer:   <xdegaye at users dot sourceforge dot net>
"
" Configure VIM to be used with pyclewn and netbeans.
"
if exists("s:did_start")
    finish
endif
let s:did_start = 1

" The following global variables define how pyclewn is started. They may be
" changed to suit your preferences.
function s:init(debugger)
    if exists("g:pyclewn_terminal")
        let s:terminal = g:pyclewn_terminal
    else
        let s:terminal = ""
    endif

    if exists("g:pyclewn_python")
        let s:python = g:pyclewn_python
    else
        let s:python = "python"
    endif

    if exists("g:pyclewn_args")
        let s:args = g:pyclewn_args
    else
        let s:args = "--window=top --maxlines=10000 --background=Cyan,Green,Magenta"
    endif

    if exists("g:pyclewn_connection")
        let s:connection = g:pyclewn_connection
    else
        if a:debugger == "gdb"
            let s:connection = "127.0.0.1:3219:changeme"
        elseif a:debugger == "pdb"
            let s:connection = "127.0.0.1:3220:changeme"
        else
            let s:connection = "127.0.0.1:3221:changeme"
        endif
    endif

    " Uncomment the following line to print full traces in a file named
    " 'logfile' for debugging purpose (or change g:pyclewn_args).
    " let s:args .= " --level=nbdebug --file=logfile"
    if s:terminal != ""
        let s:args .= " --level=info"
    endif

    let l:fixed_args = "--editor= --netbeans=" . s:connection . " --cargs="
    if s:terminal != ""
        let s:fixed = l:fixed_args
    else
        let s:fixed = "--daemon " . l:fixed_args
    endif
endfunction

" Run the 'Cinterrupt' command to open the console.
function s:interrupt(args)
    " find the prefix
    let argl = split(a:args)
    let prefix = "C"
    let idx = index(argl, "-x")
    if idx == -1
        let idx = index(argl, "--prefix")
        if idx == -1
            for item in argl
                if stridx(item, "--prefix") == 0
                    let pos = stridx(item, "=")
                    if pos != -1
                        let prefix = strpart(item, pos + 1)
                    endif
                endif
            endfor
        endif
    endif

    if idx != -1 && len(argl) > idx + 1
        let prefix = argl[idx + 1]
    endif

    " hack to prevent Vim being stuck in the command line with '--More--'
    echohl WarningMsg
    echo "About to run the 'interrupt' command."
    call inputsave()
    call input("Press the <Enter> key to continue.")
    call inputrestore()
    echohl None
    exe prefix . "interrupt"
endfunction

" Check wether pyclewn successfully wrote the script file.
function s:pyclewn_ready(filename)
    let l:cnt = 1
    echohl WarningMsg
    while l:cnt < 20
        echon "."
        let l:cnt = l:cnt + 1
        if filereadable(a:filename)
            break
        endif
        sleep 200m
    endwhile
    echohl None
    if !filereadable(a:filename)
        throw "Error: pyclewn failed to start.\n\n"
    endif
    call s:info("Creation of vim script file \"" . a:filename . "\": OK.\n")
endfunction

" Start pyclewn and vim netbeans interface.
function s:start(args)
    if !exists(":nbstart")
        call s:error("Error: the ':nbstart' vim command does not exist.")
        return
    endif
    if has("netbeans_enabled")
        call s:error("Error: netbeans is already enabled and connected.")
        return
    endif
    if !executable(s:python)
        call s:error("Error: '" . s:python . "' cannot be found or is not an executable.")
        return
    endif
    let l:tmpfile = tempname()

    " Remove the Console and the list buffers from the previous session.
    if bufexists("(clewn)_console")
        bwipeout (clewn)_console
    endif
    for l:b in ["Variables", "Breakpoints", "Backtrace", "Threads"]
        let l:bufname = "(clewn)_" . tolower(l:b)
        if bufexists(l:bufname)
            exe "bwipeout " . l:bufname
        endif
    endfor

    " Start pyclewn and netbeans.
    call s:info("Starting pyclewn.\n")
    let l:run_pyclewn = s:python . " -m clewn " . s:fixed . l:tmpfile . " " . a:args
    if s:terminal == ""
        exe "silent !" . l:run_pyclewn . " &"
    else
        let l:run_terminal = join(split(s:terminal, ","), " ")
        exe "silent !" . l:run_terminal . " sh -c '" . l:run_pyclewn . " || sleep 600' &"
    endif

    call s:info("Running nbstart, <C-C> to interrupt.\n")
    call s:pyclewn_ready(l:tmpfile)
    exe "nbstart :" . s:connection

    " source vim script
    if has("netbeans_enabled")
        " the pyclewn generated vim script is sourced only once
        if ! exists("s:source_once")
            let s:source_once = 1
            exe "source " . l:tmpfile
        endif
        call s:info("The netbeans socket is connected.\n")
        let argl = split(a:args)
        if index(argl, "pdb") == len(argl) - 1
            call s:interrupt(a:args)
        endif
    else
        throw "Error: the netbeans socket could not be connected."
    endif
endfunction

function pyclewn#start#StartClewn(...)
    let l:debugger = "gdb"
    if a:0 != 0
        let l:debugger = a:1
    endif
    call s:init(l:debugger)

    let l:args = s:args
    if a:0 != 0
        if index(["gdb", "pdb", "simple"], a:1) == -1
            call s:error("Unknown debugger '" . a:1 . "'.")
            return
        endif
        if a:0 > 1
            let l:args .= " --args \"" . join(a:000[1:], ' ') . "\""
        endif
        let l:args .= " " . a:1
    endif

    try
        call s:start(l:args)
    catch /^Vim:Interrupt$/
        return
    catch /.*/
        call s:info("The 'Pyclewn' command has been aborted.\n")
        let l:err = v:exception . "\n"
        let l:err .= "To get the cause of the problem set the global variable"
        let l:err .= " 'pyclewn_terminal' to:\n"
        let l:err .= ":let g:pyclewn_terminal = \"xterm, -e\"\n"
        call s:error(l:err)
        " The vim console screen is garbled, redraw the screen.
        if !has("gui_running")
            redraw!
        endif
        " Clear the command line.
        echo "\n"
    endtry
endfunction

function s:info(msg)
    echohl WarningMsg
    echo a:msg
    echohl None
endfunction

function s:error(msg)
    echohl ErrorMsg
    echo a:msg
    call inputsave()
    call input("Press the <Enter> key to continue.")
    call inputrestore()
    echohl None
endfunction
autoload/pyclewn/buffers.vim	[[[1
228
" vi:set ts=8 sts=4 sw=4 et tw=80:
" Pyclewn run time file.
" Maintainer:   <xdegaye at users dot sourceforge dot net>
"
" Manage pyclewn buffers.
"
if exists("s:did_buffers")
    finish
endif
let s:did_buffers = 1

"---------------------   AUTOLOAD FUNCTIONS   ---------------------

" Display one of the pyclewn buffers in a window. The function is triggered by a
" 'BufAdd' autocommand. The function is also called directly with 'name' as
" "(clewn)_console" just before the first 'C' command when the buffer list is
" empty (to workaround a problem with Vim that fails to send netbeans events
" when the buffer list is empty).
"   'name':     the pyclewn buffer name.
"   'location': the value of the '--window' option, i.e. "top", "bottom",
"               "left", "right" or "none".
function pyclewn#buffers#CreateWindow(name, location)
    if a:name == "(clewn)_empty"
        return
    endif
    if exists("*Pyclewn_CreateWindow")
        call Pyclewn_CreateWindow(a:name, a:location)
        return
    endif
    call s:create_window(a:name, a:location)
endfunction

" Display the '(clewn)_variables' buffer in a window, split if needed. The
" function is called before the 'Cdbgvar' command is executed.
function pyclewn#buffers#DbgvarSplit()
    if exists("*Pyclewn_DbgvarSplit")
        call Pyclewn_DbgvarSplit()
        return
    endif
    call s:split_clewnbuffer("(clewn)_variables", "")
endfunction

" Display the frame source code in a window. The function is called after the
" <CR> key or the mouse is used in a '(clewn)_backtrace' window. The line number
" is not available (to avoid screen blinks) in this window, but the ensuing
" 'Cframe' command will automatically move the cursor to the right place.
"   'fname': the source code full path name.
function pyclewn#buffers#GotoFrame(fname)
    if exists("*Pyclewn_GotoFrame")
        call Pyclewn_GotoFrame(a:fname)
        return
    endif
    call s:split_source(a:fname, "")
endfunction

" Display the breakpoint source code in a window. The function is called after
" the <CR> key or the mouse is used in a '(clewn)_breakpoints' window.
"   'fname': the source code full path name.
"   'lnum':  the source code line number.
function pyclewn#buffers#GotoBreakpoint(fname, lnum)
    if exists("*Pyclewn_GotoBreakpoint")
        call Pyclewn_GotoBreakpoint(a:fname, a:lnum)
        return
    endif
    call s:split_source(a:fname, a:lnum)
endfunction

"-------------------   END AUTOLOAD FUNCTIONS   -------------------

" The '(clewn)_empty' buffer is used here to workaround the problem that
" BufWinLeave auto commands are never triggered when the clewn buffer is loaded
" in a window whose current buffer is a netbeans created file.
function s:create_window(name, location)
    if a:name == "(clewn)_console"
        " When the buffer list is empty, do not split the window.
        if bufname("%") == ""
            exe "edit (clewn)_empty"
        else
            call s:split_clewnbuffer(a:name, a:location)
        endif
        return
    endif

    if a:name == "(clewn)_variables" || a:location != "top"
        return
    endif

    " Search for any existing list buffer window.
    let l:list_buffers = {'breakpoints':1, 'backtrace':2, 'threads':3}
    let l:gotit = 0
    for l:buf in keys(l:list_buffers)
        let l:name = "(clewn)_" . l:buf
        let l:nr = bufwinnr(l:name)
        if l:nr != -1
            if l:name == a:name
                return
            endif
            let l:gotit = 1
        endif
        if l:name == a:name
            let l:bufnr = l:list_buffers[l:buf]
        endif
    endfor

    let l:prevbuf_winnr = bufwinnr(bufname("%"))
    let l:count = 0
    if bufwinnr("(clewn)_console") == 1
        let l:count = 1
    endif

    if ! l:gotit
        " Create the 3 windows on the first BufAdd event of a list buffer.
        wincmd w
        if l:count
            exe (&previewheight - 4) . "split"
            wincmd w
        else
            4split
        endif
        exe "edit (clewn)_empty"
        vsplit | vsplit
    endif

    " Edit the new buffer.
    let l:bufnr = l:bufnr + l:count
    exe l:bufnr . "wincmd w"
    exe "edit " . a:name
    setlocal nowrap

    exe l:prevbuf_winnr . "wincmd w"
endfunction

" Split a window and display a buffer with previewheight.
function s:split_clewnbuffer(fname, location)
    if a:location == "none"
        return
    endif

    " The window does not exist.
    let l:nr = bufwinnr(a:fname)
    if l:nr == -1
        call s:split_location(a:fname, a:location)
    endif

    " Split the window (when the only window) this is required to prevent Vim
    " display toggling between clewn console and the last buffer where the
    " cursor was positionned (clewn does not know that this buffer is not
    " anymore displayed).
    if winnr("$") == 1
        call s:split_location("", a:location)
    endif
endfunction

" Split a window and return to the initial window,
" if 'location' is not ''
"   'location' may be: '', 'top', 'bottom', 'left' or 'right'.
function s:split_location(fname, location)
    let l:nr = 1
    let l:split = "split"
    let l:spr = &splitright
    let l:sb = &splitbelow
    set nosplitright
    set nosplitbelow
    let l:prevbuf_winnr = bufwinnr(bufname("%"))
    if winnr("$") == 1 && (a:location == "right" || a:location == "left")
        let l:split = "vsplit"
        if a:location == "right"
            set splitright
        else
            let l:prevbuf_winnr = 2
        endif
    else
        if a:location == "bottom"
            let l:nr = winnr("$")
            set splitbelow
        else
            let l:prevbuf_winnr = l:prevbuf_winnr + 1
        endif
        if a:location != ""
            exe l:nr . "wincmd w"
        endif
    endif
    let l:nr = bufnr(a:fname)
    if l:nr != -1
        exe &previewheight . l:split
        exe l:nr . "buffer"
        setlocal nowrap
    else
        exe &previewheight . l:split . " " . a:fname
        setlocal nowrap
    endif
    let &splitright = l:spr
    let &splitbelow = l:sb
    exe l:prevbuf_winnr . "wincmd w"
endfunc

function s:split_source(fname, lnum)
    let l:nr = bufwinnr(a:fname)
    if l:nr != -1
        exe l:nr . "wincmd w"
        if a:lnum != ""
            call cursor(a:lnum, 0)
        endif
        return
    endif

    " Search for a source code window.
    let l:count = winnr('$')
    let l:nr = 1
    while l:nr <= l:count
        if bufname(winbufnr(l:nr)) !~# "^(clewn)_"
            exe l:nr . "wincmd w"
            break
        endif
        let l:nr = l:nr + 1
    endwhile

    " Split the window.
    exe &previewheight . "split"
    if l:nr > l:count
        wincmd w
    endif
    exe "edit " . a:fname
    if a:lnum != ""
        call cursor(a:lnum, 0)
    endif
endfunction

doc/pyclewn.txt	[[[1
1181
*pyclewn.txt*                                   Last change: 2015 February 27


                            PYCLEWN USER MANUAL

The Pyclewn user guide                              *pyclewn*

1. Introduction                                     |pyclewn-intro|
2. Starting pyclewn                                 |:Pyclewn|
3. Options                                          |pyclewn-options|
4. Using pyclewn                                    |pyclewn-using|
5. Gdb                                              |pyclewn-gdb|
6. Pdb                                              |pyclewn-pdb|
7. Key mappings                                     |pyclewn-mappings|
8. Watched variables                                |pyclewn-variables|
9. Pyclewn windows                                  |pyclewn-windows|


==============================================================================
1. Introduction                                     *pyclewn-intro*


Pyclewn is a python program that allows the use of Vim as a front end to a
debugger. Pyclewn supports the gdb and the pdb debuggers. Pyclewn uses the
netbeans protocol to control Vim.

The debugger output is redirected to a Vim buffer. The debugger commands are
mapped to Vim user-defined commands with a common letter prefix (the default
is the |C| letter), and with Vim command completion available on the commands
and their first argument.


Pyclewn currently supports the following debuggers:

    * gdb:      version 7.1 and above, pyclewn uses the gdb MI interface.

    * pdb:      the Python debugger.

    * simple:   a fake debugger implemented in python to test pyclewn
                internals.


Pyclewn provides the following features:
---------------------------------------
* A debugger command can be mapped in Vim to a key sequence using Vim key
  mappings. This allows, for example, to set/clear a breakpoint or print a
  variable value at the current cursor or mouse position by just hitting a
  key.

* A sequence of gdb commands can be run from a Vim script when the
  |async-option| is set. This may be useful in a key mapping.

* Breakpoints and the line in the current frame are highlighted in the source
  code. Disabled breakpoints are noted with a different highlighting color.
  Pyclewn automatically finds the source file for the breakpoint if it exists,
  and tells Vim to load and display the file and highlight the line.

* The value of an expression or variable is displayed in a balloon in gvim
  when the mouse pointer is hovering over the selected expression or the
  variable.

* Similarly to gdb, one may attach to a running python process with the pdb
  debugger, interrupt the process, manage a debugging session and terminate
  the debugging session by detaching from the process. A new debugging session
  may be conducted later on this same process, possibly from another Vim
  instance.

* A gdb expression can be watched in a Vim window. The expression value is
  updated and highlighted whenever it has changed. When the expression is a
  structure or class instance, it can be expanded (resp. folded) to show
  (resp. hide) its members and their values. This feature is only available
  with gdb.

* With gdb, one can jump with the <CR> key or the mouse to the corresponding
  source code line from a "(clewn)_breakpoints" window or switch to the
  corresponding frame from a "(clewn)_backtrace" window, or switch to the
  correponding thread with a "(clewn)_threads" window.

* The |project-command| saves the current gdb settings to a project file that
  may be sourced later by the gdb "source" command. These settings are the
  working directory, the debuggee program file name, the program arguments and
  the breakpoints. The sourcing and saving of the project file can be
  automated to occur on each gdb startup and termination, whith the
  |project-file| command line option. The ``project`` command is currently only
  available with gdb.

* Vim command completion on the commands and their first argument.


The remaining sections of this manual are:
-----------------------------------------
    2. |:Pyclewn| explains how to start pyclewn.

    3. |pyclewn-options| lists pyclewn options and their usage.

    4. |pyclewn-using| explains how to use the pyclewn features common to all
        supported debuggers.

    5. |pyclewn-gdb| details the topics relevant to pyclewn and the gdb
        debugger.

    6. |pyclewn-pdb| details the topics relevant to pyclewn and the pdb
        debugger.

    7. |pyclewn-mappings| lists the pyclewn key mappings and how to use them.

    8. |pyclewn-variables| explains how to use the variable debugger window
        with gdb.

    9. |pyclewn-windows| list the windows management functions that can
        be implemented to customize the behavior of pyclewn.

==============================================================================
2. Starting pyclewn                                 *:Pyclewn*


Start pyclewn from Vim:
-----------------------
The |:Pyclewn| Vim command requires at least Vim 7.3. Start pyclewn from Vim
with the command: >

    :Pyclewn [debugger] [args]

Where "debugger" is either "gdb" or "pdb" and defaults to "gdb", and "args"
are the arguments passed to the debugger. Please note that double quotes in
|:Pyclewn| arguments are not supported. For example, to debug
"foobar foobar_arg" with gdb, run the commmand: >

    :Pyclewn gdb --args foobar foobar_arg

The gdb debugger is started upon the invocation of the first gdb command. For
example, load foobar with the gdb command "file" and start gbd by typing on
Vim command line: >

    :Cfile foobar

To just start gdb with a command that does not have any effect: >

    :Cecho

To terminate pyclewn and the Vim netbeans interface, run the following
command: >

    :nbclose

To know if the netbeans interface is connected, run the following command: >

    :echo has("netbeans_enabled")

The |:Pyclewn| command does the following:

    * spawn pyclewn
    * start the Vim netbeans interface and connect it to pyclewn
    * source a script automatically generated by pyclewn containing utility
      functions and the debugger commands as Vim commands


Global variables         *pyclewn_python* *pyclewn_args* *pyclewn_connection*
----------------
When starting pyclewn from Vim, the pyclewn command line arguments and
connection details may be set with the following global variables:

    * |pyclewn_terminal|: start pyclewn in a terminal (see next section)
    * |pyclewn_python|: the program to use as python
    * |pyclewn_args|: pyclewn arguments
    * |pyclewn_connection|: pyclewn connection parameters

When those global variables are not set, pyclewn is spawned with the
following values:

    * pyclewn_terminal: ""
    * pyclewn_python: "python"
    * gdb pyclewn_connection: "localhost:3219:changeme"
    * pdb pyclewn_connection: "localhost:3220:changeme"
    * pyclewn_args:
        "--window=top --maxlines=10000 --background=Cyan,Green,Magenta"


Trouble shooting ':Pyclewn'                          *pyclewn_terminal*
---------------------------
Set the |pyclewn_terminal| global variable to have pyclewn started in a
terminal instead of as a daemon, and visualize the log traces at the "info"
log level.|pyclewn_terminal| value is a comma separated list of the arguments
needed to start the terminal and a program running in this terminal. For
example, with xterm: >

    :let g:pyclewn_terminal = "xterm, -e"

When pyclewn fails (its exit code is not zero), the terminal remains open for
10 minutes so that the log traces may be examined. Type <Ctl-C> to close the
terminal.


Start pyclewn from a shell:
---------------------------
Start pyclewn with: >

    python -m clewn [options] [debugger]

"debugger" defaults to "gdb".

So pyclewn with the gdb debugger is simply started as: >

    python -m clewn


Start pyclewn from a Vim script:
--------------------------------
For example, to debug with gdb a program named foobar, write the following
script.vim: >

    " A function to test the ':Pyclewn' command in a script.
    function PyclewnScripting()
        " This function can be called while the previous netbeans session is
        " still open.
        if has("netbeans_enabled")
            echo "Error: netbeans is already connected."
            call input("Press the <Enter> key to continue.")
            return
        endif

        let g:pyclewn_args="--gdb=async"
        Pyclewn gdb
        Cfile foobar
        Cstart
    endfunc

and run the command: >

    gvim -S script.vim -c "call PyclewnScripting()"

==============================================================================
3. Options                                          *pyclewn-options*


The pyclewn options can be set:

    * on the command line

    * with the |pyclewn_python|, |pyclewn_args| and |pyclewn_connection| Vim
      global variables when starting pyclewn with the |:Pyclewn| command

    * as the keyword parameters of the pdb function


Options:
  --version                   show program's version number and exit
  -h, --help                  show this help message and exit
  -g PARAM_LIST, --gdb=PARAM_LIST
                              set gdb PARAM_LIST
  -d, --daemon                run as a daemon (default 'False')
  -e EDITOR, --editor=EDITOR  set Vim program to EDITOR
  -c ARGS, --cargs=ARGS       set Vim arguments to ARGS
  -p PGM, --pgm=PGM           set the debugger pathname to PGM
  -a ARGS, --args=ARGS        set the debugger arguments to ARGS
  --terminal=TERMINAL         set the terminal to use with the inferiortty
                              command (default 'xterm,-e')
  --run                       allow the debuggee to run after the pdb() call
                              (default 'False')
  --tty=TTY                   use TTY for input/output by the python script
                              being debugged (default '/dev/null')
  -w LOCATION, --window=LOCATION
                              open the debugger console window at LOCATION
                              which may be one of (top, bottom, left,
                              right, none), the default is top
  -m LNUM, --maxlines=LNUM    set the maximum number of lines of the debugger
                              console window to LNUM (default 10000 lines)
  -x PREFIX, --prefix=PREFIX  set the commands prefix to PREFIX (default 'C')
  -b COLORS, --background=COLORS
                              COLORS is a comma separated list of the three
                              colors of the breakpoint enabled, breakpoint
                              disabled and frame sign background colors, in
                              this order (default 'Cyan,Green,Magenta')
  -n CONN, --netbeans=CONN    set netBeans connection parameters to CONN with
                              CONN as 'host[:port[:passwd]]'
  -l LEVEL, --level=LEVEL     set the log level to LEVEL: critical, error,
                              warning, info, debug, nbdebug (default error)
  -f FILE, --file=FILE        set the log file name to FILE


The full description of pyclewn options follows:
------------------------------------------------
--version           Show program's version number and exit.

-h
--help              Show this help message and exit.

-g {PARAM_LIST}
--gdb={PARAM_LIST}  The PARAM_LIST option parameter is a comma separated list
                    of parameters and is mandatory when the option is present.
                    So, to run gdb with no specific parameter, the following
                    commands are equivalent: >

                        python -m clewn
                        python -m clewn -g ""
                        python -m clewn --gdb=
.
                    There are two optional parameters:

                        * the "async" keyword sets the |async-option|
                        * the project file name sets the |project-file|

                    The project file name can be an absolute pathname, a
                    relative pathname starting with '.' or a home relative
                    pathname starting with '~'. The directory of the project
                    file name must be an existing directory.
                    For example: >

                        python -m clewn --gdb=async,./project_name

-d
--daemon            Run as a daemon (default 'False'): pyclewn is detached
                    from the terminal from where it has been launched, which
                    means that this terminal cannot be used as a controlling
                    terminal for the program to debug, and cannot be used for
                    printing the pyclewn logs as well.

-e {EDITOR}
--editor={EDITOR}   Set Vim program to EDITOR. EDITOR must be in one
                    of the directories listed in the PATH environment
                    variable, or the absolute pathname of the Vim executable.
                    When this command line option is not set, pyclewn uses the
                    value of the EDITOR environment variable, and if this
                    environment variable is not set either, then pyclewn
                    defaults to using "gvim" as the name of the program to
                    spawn. Vim is not spawned by pyclewn when this option is
                    set to an empty string or when the debugger is pdb.

-c {ARGS}
--cargs={ARGS}      Set the editor program arguments to ARGS, possibly double
                    quoted (same as option --args).

-p {PGM}
--pgm={PGM}         Set the debugger program to PGM. PGM must be in one of the
                    directories listed in the PATH environment variable.

-a {ARGS}
--args={ARGS}       Set the debugger program arguments to ARGS. These
                    arguments may be double quoted. For example, start gdb
                    with the program foobar and "this is foobar argument" as
                    foobar's argument: >

                    python -m clewn -a '--args foobar "this is foobar argument"'

--terminal=TERMINAL Set the terminal to use with the inferiortty command for
                    running gdb or pdb inferior (default 'xterm,-e'). The
                    option is a comma separated list of the arguments needed
                    to start the terminal and a program running in this
                    terminal.

--run               By default the python debuggee is stopped at the first
                    statement after the call to pdb(). Enabling this option
                    allows the debuggee to run after the call to pdb().

--tty={TTY}         Use TTY for input/output by the python script being
                    debugged. The default is "/dev/null".

-w {LOCATION}
--window={LOCATION} The debugger console window pops up at LOCATION, which may
                    be one of top, bottom, left, right or none. The default is
                    top.  In the left or right case, the window pops up on the
                    left (resp. right) if there is only one window currently
                    displayed, otherwise the debugger window is opened at the
                    default top. When LOCATION is none, the automatic display
                    of the console is disabled.

-m {LNUM}
--maxlines={LNUM}   Set the maximum number of lines of the debugger console
                    window to LNUM (default 10000 lines). When the number of
                    lines in the buffer reaches LNUM, 10% of LNUM first lines
                    are deleted from the buffer.

-x {PREFIX}
--prefix={PREFIX}   Set the user defined Vim commands prefix to PREFIX
                    (default |C|). The prefix may be more than one letter
                    long. The first letter must be upper case.

-b {COLORS}
--background={COLORS}
                    COLORS is a comma separated list of the three colors of
                    the breakpoint enabled, breakpoint disabled and frame sign
                    background colors, in this order (default
                    'Cyan,Green,Magenta'). The color names are case sensitive.
                    See |highlight-ctermbg| for the list of the valid color
                    names.

                    This option has no effect when Vim version is vim72 or
                    older.

-n {CONN}
--netbeans={CONN}   Set netBeans connection parameters to CONN with CONN as
                    'host[:port[:passwd]]', (the default is
                    '127.0.0.1:3219:changeme' for gdb and
                    '127.0.0.1:3220:changeme' for pdb). Pyclewn listens on
                    host:port, with host being a name or the IP address of one
                    of the local network interfaces in standard dot notation.

-l {LEVEL}
--level={LEVEL}     Set the log level to LEVEL: critical, error, warning, info,
                    debug or nbdebug (default critical). Level nbdebug is very
                    verbose and logs all the netbeans pdu as well as all the
                    debug traces. Critical error messages are printed on
                    stderr. No logging is done on stderr (including critical
                    error messages) when the "--level" option is set to
                    something else than "critical" and the "--file" option is
                    set.

-f {FILE}
--file={FILE}       Set the log file name to FILE.

==============================================================================
4. Using pyclewn                                            *pyclewn-using*


Commands:                                           *Ccommand* *C*
---------
The prefix letter |C| is the default Vim command prefix used to map debugger
commands to Vim user-defined commands. These commands are called |Ccommand| in
this manual. The prefix can be changed with a command line option.

A debugger command can be entered on Vim command line with the |C| prefix. It is
also possible to enter the command as the first argument of the |C| command. In
the following example with gdb, both methods are equivalent: >

    :Cfile /path/to/foobar
    :C file /path/to/foobar

The first method provides completion on the file name while the second one
does not.

The second method is useful when the command is a user defined command in the
debugger (user defined commands built by <define> in gdb), and therefore not a
Vim command. It is also needed for gdb command names that cannot be mapped to
a Vim command because Vim does not accept non alphanumeric characters within
command names (for example <core-file> in gdb).

To get help on the pyclewn commands, use Chelp.

Pyclewn commands can be mapped to keys, or called within a Vim script or a
menu.

Note:
The gdb debugger cannot handle requests asynchronously, so the
|async-option| must be set, when mapping a key to a sequence of commands.
With this option set, one can build for example the following mapping: >

    :map <F8> :Cfile /path/to/foobar <Bar> Cbreak main <Bar> Crun <CR>

Note:
Quotes and backslashes must be escaped on Vim command line. For example, to
print foo with a string literal as first argument to the foo function: >

    :Cprint foo(\"foobar\", 1)

And to do the same thing with the string including a new line: >

    :Cprint foo(\"foobar\\n\", 1)


Completion:
-----------
Command line completion in Vim is usually done using the <Tab> key (set by the
'wildchar' option). To get the list of all valid completion matches, type
CTRL-D. For example, to list all the debugger commands (assuming the
default |C| prefix is being used): >

    :C<C-D>

See also the 'wildmenu' option. With this option, the possible matches are
shown just above the command line and can be selected with the arrow keys.

The first argument completion of a |Ccommand| may be done on a file name or on a
list. For example with gdb, the following command lists all the gdb help
items: >

    :Chelp <C-D>

The first argument completion of the |C| command is the list of all the debugger
commands. For example, to list all the debugger commands (note the space after
the |C|): >

    :C <C-D>


Command line search:
--------------------
Use the powerful command line search capabilities of the Vim command line.
For example, you want to type again, possibly after a little editing, one of
the commands previously entered: >

    :Cprint (*(foo*)0x0123ABCD)->next->next->part1->something_else.aaa

You can get rapidly to this command by using the Vim command line window
|cmdline-window|: >

    :<CTRL-F>
    /something_else
    <CR>

or from normal mode >
    q:
    /something_else
    <CR>


Vim in a terminal
-----------------
The debuggee output is redirected to '/dev/null' when the name of the program
is "vim" or "vi". Use the |inferiortty| command to redirect the debuggee
output to a terminal.

Do not use the "--daemon" command line option when running Vim in a console.


Balloon:
--------
A variable is evaluated by the debugger and displayed in a balloon in gvim,
when the mouse pointer is hovering over the the variable. To get the
evaluation of an expression, first select the expression in visual mode in the
source code and point the mouse over the selected expression. Disable this
feature with the |Cballooneval| command.

==============================================================================
5. Gdb                                              *pyclewn-gdb*


When gdb is started, it automatically executes commands from its init file,
normally called '.gdbinit'. See the gdb documentation.


                                                    *clewn-list-buffers*

The clewn list buffers:
-----------------------
Three buffers are updated by pyclewn whenever the state of the debuggee
changes. The buffers list the breakpoints, the backtrace and the threads. From
a window loaded with one of these buffers and with the <CR> key or the mouse,
one may:

    * "(clewn)_breakpoints": jump to the corresponding source code line.

    * "(clewn)_backtrace": switch to the corresponding frame.

    * "(clewn)_threads": switch to the correponding thread.

When the "--window" option is "top" (the default), pyclewn creates three
windows at the first user command that creates a breakpoint or start the
debuggee. The windows are then populated with their corresponding list buffers
as soon as these buffers become available. See |pyclewn-windows| for a way to
disable or customize this feature.


                                                    *inferior_tty*

Debuggee standard input and output:
-----------------------------------
When starting pyclewn from a terminal and using gvim, pyclewn creates a pseudo
terminal that is the the controlling terminal of the program to debug.
Programs debugged by gdb, including those based on curses and termios such as
Vim, run in this terminal. A <Ctl-C> typed in the terminal interrupts the
debuggee.

When pyclewn is started from Vim with the |:Pyclewn| command, there is no
terminal associated with pyclewn. The |inferiortty| command provides the same
functionality as above and spawns the controlling terminal (using the
--terminal option, default xterm) of the debuggee and sets accordingly gdb
'inferior-tty' variable and the TERM environment variable. The gdb
'inferior-tty' variable MUST be set BEFORE the inferior is started.


                                                    *async-option*
Async option:
-------------
The gdb event loop is not asynchronous in most configurations, which means
that gdb cannot handle a command while the previous one is being processed and
discards it.
When gdb is run with the |async-option| set, pyclewn queues the commands in a
fifo and send a command to gdb, only when gdb is ready to process the command.
This allows the key mappings of a sequence of gdb commands. To set the
|async-option| , see |pyclewn-options|.


                                                    *gdb-keys*
List of the gdb default key mappings:
-------------------------------------
These keys are mapped after the |Cmapkeys| Vim command is run.

        CTRL-Z  send an interrupt to gdb and the program it is running
        B       info breakpoints
        L       info locals
        A       info args
        S       step
        CTRL-N  next: next source line, skipping all function calls
        F       finish
        R       run
        Q       quit
        C       continue
        W       where
        X       foldvar
        CTRL-U  up: go up one frame
        CTRL-D  down: go down one frame

cursor position: ~
        CTRL-B  set a breakpoint on the line where the cursor is located
        CTRL-K  clear all breakpoints on the line where the cursor is located

mouse pointer position: ~
        CTRL-P  print the value of the variable defined by the mouse pointer
                position
        CTRL-X  print the value that is referenced by the address whose
                value is that of the variable defined by the mouse pointer
                position


                                                    *$cdir*
Source path:
-----------
Pyclewn automatically locates the source file with the help of gdb, by using
the debugging information stored in the file that is being debugged when the
gdb "directories" setting has a component named "$cdir" (the default). This is
useful when the program to debug is the result of multiple compilation units
located in different directories.


                                                    *Csymcompletion*
Symbols completion:
-------------------
The gdb <break> and <clear> commands are set initially with file name
completion. This can be changed to completion matching the symbols of the
program being debugged, after running the |Csymcompletion| command. This is a
pyclewn command.

To minimize the number of loaded symbols and to avoid fetching the shared
libraries symbols, run the Csymcompletion command after the file is loaded
with the gdb <file> command, and before the program is run.

Note: The <break> and <clear> filename completion is not the same as gdb file
name completion for these two commands. Gdb uses the symbols found in the
program file to debug, while pyclewn uses only the file system.


                                                    *gdb-balloon*
Balloon evaluation:
-------------------
To dereference a pointer "ptr" to a structure and show its members, select an
expression that dereferences the pointer, for example "*ptr", and hover above
the selected area with the mouse. Another solution is to hover with the mouse
above "ptr" and type <C-X>.


                                      *Cproject* *project-command* *project-file*
Project file:
-------------
The pyclewn |project-command| name is "project". This command saves the current
gdb settings to a project file that may be sourced later by the gdb "source"
command.

These settings are:
    * current working directory
    * debuggee program file name
    * program arguments
    * all the breakpoints (at most one breakpoint per source line is saved)

The argument of the |project-command| is the pathname of the project file.
For example: >

    :Cproject /path/to/project

When the "--gdb" option is set with a project filename (see |pyclewn-options|),
the project file is automatically sourced when a a gdb session is started, and
the project file is automatically saved when the gdb session or Vim session,
is terminated.

Note: When gdb sources the project file and cannot set a breakpoint because,
for example, it was set in a shared library that was loaded at the time the
project file was saved, gdb ignores silently the breakpoint (see gdb help on
"set breakpoint pending").


Limitations:
------------
When setting breakpoints on an overloaded method, pyclewn bypasses the gdb
prompt for the multiple choice and sets automatically all breakpoints.

In order to set a pending breakpoint (for example in a shared library that has
not yet been loaded by gdb), you must explicitly set the breakpoint pending
mode to "on", with the command: >

    :Cset breakpoint pending on

After a "detach" gdb command, the frame sign remains highlighted because
gdb/mi considers the frame as still valid.

When answering "Abort" to a dialog after pyclewn attempts to edit a buffer and
set a breakpoint in a file already opened within another Vim session, the
breakpoint is set in gdb, but not highlighted in the corresponding buffer.
However, it is possible to |bwipeout| a buffer at any time, and load it again in
order to restore the correct highlighting of all the breakpoints in the
buffer.


Pyclewn commands:
-----------------
The |Ccommand| list includes all the gdb commands and some pyclewn specific
commands that are listed here:

                                                            *Cballooneval*
    * Cballooneval   enable or disable showing text in the Vim balloon; this
                     can also accomplished with the 'ballooneval' Vim option
                     but using |Cballooneval| allows also the mappings such as
                     CTRL-P to be still functioning after showing text in the
                     balloon has been disabled.

    *|Cdbgvar|       add a watched variable or expression to the
                     (clewn)_variables buffer.

    *|Cdelvar|       delete a watched variable from the (clewn)_variables
                     buffer.

    * Cdumprepr      print on the console pyclewn internal structures that
                     may be used for debugging pyclewn.

    *|Cfoldvar|      collapse/expand the members of a watched structure or
                     class instance.

    * Chelp          print on the console, help on the pyclewn specific
                     commands (those on this list) in addition to the help on
                     the debugger commands.
                                                    *inferiortty*
    * Cinferiortty   spawn the controlling terminal (default xterm) of the
                     debuggee and sets accordingly gdb 'inferior-tty' variable
                     and the TERM environment variable; this command  MUST be
                     issued BEFORE starting the inferior.

    * Cloglevel      print or set the log level dynamically from inside Vim.

    *|Cmapkeys|      map pyclewn keys.

    *|Cproject|      save the current gdb settings to a project file.

    *|Csetfmtvar|    set the output format of the value of a watched variable.

    * Csigint        send a <C-C> character to the debugger to interrupt the
                     running program that is being debugged; only with gdb,
                     and when pyclewn and gdb communicate over a pseudo
                     terminal.

    *|Csymcompletion|populate the break and clear commands with symbols
                     completion (only with gdb).

    * Cunmapkeys     unmap the pyclewn keys, this Vim command does not invoke
                     pyclewn.


List of illegal gdb commands:
-----------------------------
The following gdb commands cannot be run from pyclewn:

        complete
        edit
        end
        set annotate
        set confirm
        set height
        set width
        shell

==============================================================================
6. Pdb                                              *pyclewn-pdb*


Start a python script from Vim and debug it, or attach to a running python
process and start the debugging session.


Start a python script from Vim and debug it:
-------------------------------------------
To debug a python script named "script.py", run the Vim command (arg1, arg2,
... being the script.py command line arguments): >

    :Pyclewn pdb script.py arg1 arg2 ...

Or, more conveniently, debug the python script being edited in Vim as the
current buffer with: >

    :Pyclewn pdb %:p

The script is started without a controlling terminal unless the "tty" option
has been set (see below). The |Cinferiortty| command spawns a controlling
terminal (using the --terminal option that defaults to xterm) connected to a
pseudo tty, and redirects all three standard streams of the script to this
pseudo tty.

One may also redirect the script output to another tty, using the "tty" option
and setting the |pyclewn_args| Vim global variable before starting the script.
For example: >

    :let g:pyclewn_args="--tty=/dev/pts/4"

The ":Cquit" command and the Vim ":quitall" command terminate the debugging
session and the script being debugged. Both commands MUST be issued at the pdb
prompt.


Attach to a python process and debug it: >
----------------------------------------
To debug a python process after having attached to it, first insert the
following statement in the debuggee source code before starting it: >

    import clewn.vim as vim; vim.pdb()

By default, the debuggee is stopped at the first statement following the call
to vim.pdb(). To let the debuggee run instead, then use the "run" option: >

    import clewn.vim as vim; vim.pdb(run=True)

Next, attach to the process and start the debugging session by running the Vim
command: >

    :Pyclewn pdb

Notes:
The current debugging session may be terminated with the ":Cdetach" or the Vim
":quitall" command. Another debugging session with the same process can be
started later with the ":Pyclewn pdb" command.

The ":Cdetach", ":Cquit" or the Vim ":quitall" commands do not terminate the
debuggee. To kill the debuggee, issue the following command at the pdb prompt:
>
    :C import sys; sys.exit(1)

When the python process is not attached, typing two <Ctl-C> instead of one, is
needed to kill the process. This is actually a feature that allows the process
to run without any tracing overhead (before the first <Ctl-C>) when it is not
attached and no breakpoints are being set (there is still the small overhead
of the context switches between the idle clewn thread and the target thread).


Pdb commands:
-------------
The commands "interrupt", "detach" and "threadstack" are new pdb commands and
are the only commands that are available at the "[running...]" prompt when the
debuggee is running. Use the "help" command (and completion on the first
argument of the help command) to get help on each command.

The following list describes the pdb commands that are new or behave
differently from the pdb commands of the Python standard library:

                                                    *Cinterrupt*
interrupt
    This command interrupts the debuggee and is available from the
    "[running...]" prompt.

                                                    *Cinferiortty*
inferiortty
    Without argument, the pdb command "inferiortty" spawns a terminal
    connected to a pseudo tty and redirects all three standard streams to this
    pseudo tty.
    With the name of an existing pseudo tty as an argument, "inferiortty'
    redirects all three standard streams to this pseudo tty (convenient for
    re-using the same pseudo tty across multiple debugging sessions).
    This command can be issued after the script has been started.

                                                    *Cdetach*
detach
    This command terminates the debugging session by closing the netbeans
    socket. The debuggee is free to run and does not stop at the breakpoints.
    To start another debugging session, run the command: >

        :Pyclewn pdb

.   The breakpoints becomes effective again when the new session starts up.
    Available from the "[running...]" prompt and from the pdb prompt.

                                                    *Cquit*
quit
    This command terminates the debugging session by closing the netbeans
    socket, and removes the python trace function. The pyclewn thread in
    charge of handling netbeans connection terminates and it is not possible
    anymore to attach to the process. Since there is no trace function, the
    breakpoints are ineffective and the process performance is not impaired
    anymore by the debugging overhead.

    When the script has been started from Vim, this command terminates the
    script.

                                                    *Cthreadstack*
threadstack
    The command uses the sys._current_frames() function from the standard
    library to print a stack of the frames for all the threads.
    The function sys._current_frames() is available since python 2.5.
    Available from the "[running...]" prompt and from the pdb prompt.

                                                    *Cclear*
clear
    This command is the same as the Python standard library "clear" command,
    except it requires at least one parameter and therefore, it is not
    possible to clear all the breakpoints in one shot with the "clear" command
    without parameters.

the prefix alone:
    There is no "!" pdb command as in the Python standard library since Vim
    does not allow this character in a command name. However, the equivalent
    way to execute a python statement in the context of the current frame is
    with the command prefix alone, for example: >

        :C global list_options; list_options = ['-l']
        :C import sys; sys.exit(1)

.   The first word of the statement must not be a pdb command and will be
    expanded if it is an alias.

not implemented:
    The following pdb commands are not implemented: list, ll, whatis, source,
    display, undisplay, interact, run, restart.


                                                    *pdb-pdbrc*
The initialisation file .pdbrc:
-------------------------------
This file is read at initialisation and its commands are executed on startup.
See the pdb python documentation for the location of this file. Breakpoints
can be set through this file, or aliases may be defined. One useful alias
entered in the file would be for example: >

    alias kill import sys; sys.exit(1)

So that the debuggee may be killed with the command: >

    :C kill


                                                    *pdb-keys*
List of the pdb default key mappings:
-------------------------------------
These keys are mapped after the |Cmapkeys| Vim command is run.

        CTRL-Z  interrupt the pdb process
        B       list all breaks, including for each breakpoint, the number of
                times that breakpoint has been hit, the current ignore count,
                and the associated condition if any
        A       print the argument list of the current function
        S       step
        CTRL-N  next: next source line, skipping all function calls
        R       continue execution until the current function returns
        C       continue
        W       where
        CTRL-U  up: go up one frame
        CTRL-D  down: go down one frame

cursor position: ~
        CTRL-B  set a breakpoint on the line where the cursor is located
        CTRL-K  clear all breakpoints on the line where the cursor is located

mouse pointer position: ~
        CTRL-P  print the value of the selected expression defined by the
                mouse pointer position


Pyclewn commands:
-----------------
The |Ccommand| list includes pdb commands and some pyclewn specific commands
that are listed here:

    *|Cballooneval|  enable or disable showing text in the Vim balloon

    * Cdumprepr      print on the console pyclewn internal structures that
                     may be used for debugging pyclewn

    * Cloglevel      print or set the log level dynamically from inside Vim

    *|Cmapkeys|      map pyclewn keys

    * Cunmapkeys     unmap the pyclewn keys, this Vim command does not invoke
                     pyclewn


Troubleshooting:
----------------
* Pyclewn error messages can be logged in a file with the "--file" option.
  When starting the debuggee from Vim, use the |pyclewn_args| Vim global
  variable before starting the script: >

    :let g:pyclewn_args="--file=/path/to/logfile"

When attaching to a python process, use the corresponding keyword argument: >

    import clewn.vim as vim; vim.pdb(file='/path/to/logfile')


* To conduct two debugging sessions simultaneously (for example when debugging
  pyclewn with pyclewn), change the netbeans socket port with the
  |pyclewn_connection| Vim global variable before starting the script: >

    :let g:pyclewn_connection="localhost:3222:foo"

And change the corresponding keyword argument: >

    import clewn.vim as vim; vim.pdb(netbeans='localhost:3222:foo')


Limitations:
------------
The |Cinterrupt| command does not properly interrupt the input() Python
function. Workaround: after the Cinterrupt command has been issued while at
the prompt displayed by input(), enter some data to allow the input() function
to complete execution of its C code implementation, this allows pdb to gain
control when back in python code and to stop.

==============================================================================
7. Key mappings                                     *pyclewn-mappings*


All |Ccommand| can be mapped to Vim keys using the Vim |:map-commands|.
For example, to set a breakpoint at the current cursor position: >

    :map <F8> :exe "Cbreak " . expand("%:p") . ":" . line(".")<CR>

Or to print the value of the variable under the cursor: >

    :map <F8> :exe "Cprint " . expand("<cword>") <CR>


                                                    *Cmapkeys*
Pyclewn key mappings:
---------------------
This section describes another mapping mechanism where pyclewn maps Vim keys
by reading a configuration file. This is done when the |Cmapkeys| Vim command
is run. The pyclewn keys mapping is mostly useful for the pyclewn casual user.
When the configuration file cannot be found, pyclewn sets the default key
mappings. See |gdb-keys| for the list of default key mappings for gdb and
|pdb-keys| for the list of default key mappings for pdb.

Please note that pyclewn relies on the Vim |balloon-eval| feature to get the
text under the mouse position when expanding the ${text} macro. This feature
is not available with Vim console. So in this case you must build your own
key mapping as in the above example.

The configuration file is named .pyclewn_keys.{debugger}, where debugger is
the name of the debugger. The default placement for this file is
$CLEWNDIR/.pyclewn_keys.{debugger}, or $HOME/.pyclewn_keys.{debugger}.

To customize pyclewn key mappings copy the configurations files found in the
distribution to the proper directory: >

    cp runtime/.pyclewn_keys.gdb        $CLEWNDIR

or >

    cp runtime/.pyclewn_keys.gdb        $HOME

The comments in the configuration file explain how to customize the key
mappings.

Copy these files to the $CLEWNDIR or $HOME directory, and customize the key
mappings.

==============================================================================
8. Watched variables                                *pyclewn-variables*


The Watched Variables feature is available with the gdb debugger. The Vim
watched variables buffer is named "(clewn)_variables".

                                                    *Cdbgvar*
The |Cdbgvar| command is used to create a gdb watched variable in the Variables
buffer from any valid expression. A valid expression is an expression that is
valid in the current frame. The argument of the |Cdbgvar| pyclewn command is
the expression to be watched. For example, to create a watched variable for
the expression "len - max":
>
    :Cdbgvar len - max

Upon creation, the watched variable is given a name by gdb, for example: >
    <var1>
The watched variables buffer, "(clewn)_variables", is created upon creation of
the first watched variable.

                                                    *Cfoldvar*
When the watched variable is a structure or class instance, it can be expanded
and collapsed with the |Cfoldvar| command to display all its members and their
values as children watched variables. The argument of the |Cfoldvar| command
is the line number of the watched variable to expand, in the watched variable
window. For example: >

    :Cfoldvar 1

The |Cfoldvar| command is mapped to <CR> and to the double-click on the mouse
left-button so that it is easy to expand/collapse the tree with the mouse or
<CR> key.


                                                    *Cdelvar*
A gdb watched variable can be deleted with the |Cdelvar| pyclewn command.
The argument of the |Cdelvar| command is the name of the variable as given by
gdb upon creation.
For example: >

    :Cdelvar var1

When the watched variable is a structure or class instance and it has been
expanded, all its children are also deleted.


                                                    *Csetfmtvar*
Set the output format of the value of the watched variable <name>
to be <format>: >

    :Csetfmtvar <name> <format>

Parameter <name> is the gdb/mi name of the watched variable or one of its
children.
Parameter <format> is one of the strings in the following list:

    {binary | decimal | hexadecimal | octal | natural}

The "natural" format is the default format chosen automatically based on the
variable type (like "decimal" for an int, "hexadecimal" for pointers, etc.).
For a variable with children, the format is set only on the variable itself,
and the children are not affected.

Note: The setting of the format of a child watched variable is lost after
folding one of its parents (because the child is actually not watched anymore
by gdb after the folding).


Highlighting:
-------------
When the value of a watched variable has changed, it is highlighted with the
"Special" highlight group.

When a watched variable becomes out of scope, it is highlighted with the
"Comment" highlight group.

The foreground and background colors used by these highlight groups are setup
by the |:colorscheme| currently in use.

==============================================================================
9. Pyclewn windows                                *pyclewn-windows*


All the Vim functions defined by pyclewn that manage windows can be overriden.
When the following functions are defined, they will be called by pyclewn
instead of calling the corresponding functions defined in
"autoload/pyclewn/buffers.vim": >


    Pyclewn_CreateWindow(name, location)
        Display one of the pyclewn buffers in a window. The function is
        triggered by a 'BufAdd' autocommand. The function is also called
        directly with 'name' as "(clewn)_console" just before the first 'C'
        command when the buffer list is empty (to workaround a problem with
        Vim that fails to send netbeans events when the buffer list is empty).

          'name':     the pyclewn buffer name.
          'location': the value of the '--window' option, i.e. "top", "bottom",
                      "left", "right" or "none".

    Pyclewn_DbgvarSplit()
        Display the '(clewn)_variables' buffer in a window, split if needed.
        The function is called before the 'Cdbgvar' command is executed.

    Pyclewn_GotoFrame(fname)
        Display the frame source code in a window. The function is called
        after the <CR> key or the mouse is used in a '(clewn)_backtrace'
        window. The line number is not available (to avoid screen blinks) in
        this window, but the ensuing 'Cframe' command will automatically move
        the cursor to the right place.

          'fname': the source code full path name.

    Pyclewn_GotoBreakpoint(fname, lnum)
        Display the breakpoint source code in a window. The function is called
        after the <CR> key or the mouse is used in a '(clewn)_breakpoints'
        window.

          'fname': the source code full path name.
          'lnum':  the source code line number.

==============================================================================
vim:tw=78:ts=8:ft=help:norl:et:
plugin/pyclewn.vim	[[[1
11
" Pyclewn run time file.
" Maintainer:   <xdegaye at users dot sourceforge dot net>

" Enable balloon_eval.
if has("balloon_eval")
    set ballooneval
    set balloondelay=100
endif

" The 'Pyclewn' command starts pyclewn and vim netbeans interface.
command -nargs=* -complete=file Pyclewn call pyclewn#start#StartClewn(<f-args>)
syntax/clewn_variables.vim	[[[1
25
" Vim syntax file
" Language:	debugger variables window syntax file
" Maintainer:	<xdegaye at users dot sourceforge dot net>
" Last Change:	Oct 8 2007

if exists("b:current_syntax")
    finish
endif

syn region dbgVarChged display contained matchgroup=dbgIgnore start="={\*}"ms=s+1 end="$"
syn region dbgDeScoped display contained matchgroup=dbgIgnore start="={-}"ms=s+1 end="$"
syn region dbgVarUnChged display contained matchgroup=dbgIgnore start="={=}"ms=s+1 end="$"

syn match dbgItem display transparent "^.*$"
    \ contains=dbgVarUnChged,dbgDeScoped,dbgVarChged,dbgVarNum

syn match dbgVarNum display contained "^\s*\d\+:"he=e-1

high def link dbgVarChged   Special
high def link dbgDeScoped   Comment
high def link dbgVarNum	    Identifier
high def link dbgIgnore	    Ignore

let b:current_syntax = "clewn_variables"

macros/.pyclewn_keys.gdb	[[[1
49
# .pyclewn_keys.gdb file
#
# The default placement for this file is $CLEWNDIR/.pyclewn_keys.gdb, or
# $HOME/.pyclewn_keys.gdb
#
# Key definitions are of the form `KEY:COMMAND'
# where the following macros are expanded:
#    ${text}:   the word or selection below the mouse
#    ${fname}:  the current buffer full pathname
#    ${lnum}:   the line number at the cursor position
#
# All characters following `#' up to the next new line are ignored.
# Leading blanks on each line are ignored. Empty lines are ignored.
#
# To tune the settings in this file, you will have to uncomment them,
# as well as change them, as the values on the commented-out lines
# are the default values. You can also add new entries. To remove a
# default mapping, use an empty GDB command.
#
# Supported key names:
#       . key function: F1 to F20
#             e.g., `F11:continue'
#       . modifier (C-,S-,M-) + function key
#             e.g., `C-F5:run'
#       . modifier (or modifiers) + character
#             e.g., `S-Q:quit', `C-S-B:info breakpoints'
#
# Note that a modifier is required for non-function keys. So it is not possible
# to map a lower case character with this method (use the Vim 'map' command
# instead).
#
# C-B : break "${fname}":${lnum} # set breakpoint at current line
# C-D : down
# C-K : clear "${fname}":${lnum} # clear breakpoint at current line
# C-N : next
# C-P : print ${text}            # print value of selection at mouse position
# C-U : up
# C-X : print *${text}           # print value referenced by word at mouse position
# C-Z : sigint                   # kill the inferior running program
# S-A : info args
# S-B : info breakpoints
# S-C : continue
# S-F : finish
# S-L : info locals
# S-Q : quit
# S-R : run
# S-S : step
# S-W : where
# S-X : foldvar ${lnum}          # expand/collapse a watched variable
macros/.pyclewn_keys.pdb	[[[1
44
# .pyclewn_keys.pdb file
#
# The default placement for this file is $CLEWNDIR/.pyclewn_keys.pdb, or
# $HOME/.pyclewn_keys.pdb
#
# Key definitions are of the form `KEY:COMMAND'
# where the following macros are expanded:
#    ${text}:   the word or selection below the mouse
#    ${fname}:  the current buffer full pathname
#    ${lnum}:   the line number at the cursor position
#
# All characters following `#' up to the next new line are ignored.
# Leading blanks on each line are ignored. Empty lines are ignored.
#
# To tune the settings in this file, you will have to uncomment them,
# as well as change them, as the values on the commented-out lines
# are the default values. You can also add new entries. To remove a
# default mapping, use an empty GDB command.
#
# Supported key names:
#       . key function: F1 to F20
#             e.g., `F11:continue'
#       . modifier (C-,S-,M-) + function key
#             e.g., `C-F5:run'
#       . modifier (or modifiers) + character
#             e.g., `S-Q:quit', `C-S-B:info breakpoints'
#
# Note that a modifier is required for non-function keys. So it is not possible
# to map a lower case character with this method (use the Vim 'map' command
# instead).
#
# C-B : break "${fname}:${lnum}" # set breakpoint at current line
# C-D : down
# C-K : clear "${fname}:${lnum}" # clear breakpoint at current line
# C-N : next
# C-P : p ${text}                # print value of selection at mouse position
# C-U : up
# C-Z : interrupt
# S-A : args
# S-B : break
# S-C : continue
# S-R : return
# S-S : step
# S-W : where
macros/.pyclewn_keys.simple	[[[1
38
# .pyclewn_keys.simple file
#
# The default placement for this file is $CLEWNDIR/.pyclewn_keys.simple, or
# $HOME/.pyclewn_keys.simple
#
# Key definitions are of the form `KEY:COMMAND'
# where the following macros are expanded:
#    ${text}:   the word or selection below the mouse
#    ${fname}:  the current buffer full pathname
#    ${lnum}:   the line number at the cursor position
#
# All characters following `#' up to the next new line are ignored.
# Leading blanks on each line are ignored. Empty lines are ignored.
#
# To tune the settings in this file, you will have to uncomment them,
# as well as change them, as the values on the commented-out lines
# are the default values. You can also add new entries. To remove a
# default mapping, use an empty GDB command.
#
# Supported key names:
#       . key function: F1 to F20
#             e.g., `F11:continue'
#       . modifier (C-,S-,M-) + function key
#             e.g., `C-F5:run'
#       . modifier (or modifiers) + character
#             e.g., `S-Q:quit', `C-S-B:info breakpoints'
#
# Note that a modifier is required for non-function keys. So it is not possible
# to map a lower case character with this method (use the Vim 'map' command
# instead).
#
# C-B : break ${fname}:${lnum}   # set breakpoint at current line
# C-K : clear ${fname}:${lnum}   # clear breakpoint at current line
# C-P : print ${text}            # print value of selection at mouse position
# C-Z : interrupt                # interrupt the execution of the target
# S-C : continue
# S-Q : quit
# S-S : step
