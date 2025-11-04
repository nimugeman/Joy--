// joy-- interpreter in Zig
// Zig version: 0.11.0 or later

const std = @import("std");
const assert = std.debug.assert;
const stdout = std.io.getStdOut().writer();
const stderr = std.io.getStdErr().writer();
const Allocator = std.mem.Allocator;

pub const Word = u64;
const Word_Size = @sizeOf(Word);

// Type tags based on address alignment
const TAG_INT: Word = 0b1;
const TAG_LIST: Word = 0b0000;
const TAG_SYM: Word = 0b1000;
const FLAG_LINEAR: Word = 0b0010;
const MASK_ADDR: Word = ~@as(Word, 0b1111);

const T: Word = 3;
const F: Word = 1;

// Type predicates
inline fn is_int(w: Word) bool {
    return (w & 0b1) == TAG_INT;
}

inline fn is_list(w: Word) bool {
    return w != 0 and (w & 0b1111) == TAG_LIST or (w & 0b1111) == (TAG_LIST | FLAG_LINEAR);
}

inline fn is_symbol(w: Word) bool {
    return (w & 0b1111) == TAG_SYM or (w & 0b1111) == (TAG_SYM | FLAG_LINEAR);
}

inline fn is_linear(w: Word) bool {
    return !is_int(w) and w != 0 and (w & FLAG_LINEAR) != 0;
}

inline fn get_ptr(w: Word) *[2]Word {
    return @ptrFromInt(w & MASK_ADDR);
}

inline fn get_int(w: Word) i64 {
    return @as(i64, @bitCast(w)) >> 1;
}

inline fn make_int(v: i64) Word {
    return (@as(Word, @bitCast(v)) << 1) | TAG_INT;
}

// Memory pool management
var gpa: Allocator = undefined;
const POOL_SIZE: usize = 4096;
var mem_pool_head: Word = 0;
var all_pools: std.ArrayList([POOL_SIZE][2]Word) = undefined;

fn init_mem_pool() !void {
    all_pools = std.ArrayList([POOL_SIZE][2]Word).init(gpa);
    try expand_mem_pool();
}

fn deinit_mem_pool() void {
    all_pools.deinit();
}

fn expand_mem_pool() !void {
    const new_pool = try all_pools.addOne();

    var i: usize = 0;
    while (i < POOL_SIZE - 1) : (i += 1) {
        const next_addr = @intFromPtr(&new_pool[i + 1]);
        new_pool[i][0] = next_addr;
    }
    new_pool[POOL_SIZE - 1][0] = mem_pool_head;
    mem_pool_head = @intFromPtr(&new_pool[0]);
}

fn alloc_obj() !*[2]Word {
    if (mem_pool_head == 0) {
        try expand_mem_pool();
    }
    const obj: *[2]Word = @ptrFromInt(mem_pool_head);
    mem_pool_head = obj[0];
    return obj;
}

inline fn free_obj(obj: *[2]Word) void {
    obj[0] = mem_pool_head;
    mem_pool_head = @intFromPtr(obj);
}

inline fn discard(w: Word) void {
    if (!is_linear(w)) return;
    discard_linear(w);
}

fn discard_linear(w: Word) void {
    const obj = get_ptr(w);

    if (is_list(w)) {
        discard(obj[0]);
        discard(obj[1]);
    } else if (is_symbol(w)) {
        discard(obj[1]);
    }

    free_obj(obj);
}

fn deep_copy(w: Word) !Word {
    if (is_int(w)) {
        return w;
    } else if (is_list(w)) {
        const obj = get_ptr(w);
        const new_car = try deep_copy(obj[0]);
        errdefer discard(new_car);
        const new_cdr = try deep_copy(obj[1]);
        errdefer discard(new_cdr);
        return make_list(new_car, new_cdr, true);
    } else if (is_symbol(w)) {
        const obj = get_ptr(w);
        const name = obj[0];
        const op_copy = try deep_copy(obj[1]);
        errdefer discard(op_copy);
        return make_symbol(name, op_copy, true);
    }
    return w;
}

// Constructors
const NIL: Word = 0;

fn make_list(car_val: Word, cdr_val: Word, linear: bool) !Word {
    const obj = try alloc_obj();
    obj[0] = car_val;
    obj[1] = cdr_val;
    const addr = @intFromPtr(obj);
    return if (linear) addr | FLAG_LINEAR else addr;
}

fn make_symbol(name: Word, op: Word, linear: bool) !Word {
    const obj = try alloc_obj();
    obj[0] = name;
    obj[1] = op;
    const addr = @intFromPtr(obj) + Word_Size;
    return if (linear) addr | TAG_SYM | FLAG_LINEAR else addr | TAG_SYM;
}

inline fn car(w: Word) Word {
    return get_ptr(w)[0];
}

inline fn cdr(w: Word) Word {
    return get_ptr(w)[1];
}

inline fn symbol_op(w: Word) Word {
    return get_ptr(w)[1];
}

// Symbol table and string interning
var sym_table: std.StringHashMap(Word) = undefined;
var interned_strings: std.StringHashMap(Word) = undefined;

fn intern_string(s: []const u8) !Word {
    if (interned_strings.get(s)) |w| return w;
    const copy = try gpa.dupeZ(u8, s);
    const w = @intFromPtr(copy.ptr);
    try interned_strings.put(copy[0..s.len], w);
    return w;
}

fn lookup_symbol(name: []const u8) ?Word {
    return sym_table.get(name);
}

// Stack
const STACK_REAL_MAX: usize = 1024;
const STACK_MAX: usize = 1014;
var stack: [STACK_REAL_MAX]Word = undefined;
var sp: usize = 0;
var btm: usize = 0;

fn push_unchecked(w: Word) void {
    stack[sp] = w;
    sp += 1;
}

fn push(w: Word) !void {
    stack[sp] = w;
    sp += 1;
}

fn pop() Word {
    sp -= 1;
    return stack[sp];
}

fn cleanup_stack() void {
    while (sp > 0) {
        sp -= 1;
        discard(stack[sp]);
    }
}

// Primitives
const PRIM_DUP: u64 = 1;
const PRIM_SWAP: u64 = 3;
const PRIM_POP: u64 = 5;
const PRIM_DIP: u64 = 7;
const PRIM_ROT: u64 = 9;
const PRIM_CONS: u64 = 11;
const PRIM_UNCONS: u64 = 13;
const PRIM_CAT: u64 = 15;
const PRIM_I: u64 = 17;
const PRIM_FIRST: u64 = 19;
const PRIM_REST: u64 = 21;
const PRIM_IFTE: u64 = 23;
const PRIM_BRANCH: u64 = 25;
const PRIM_ADD: u64 = 27;
const PRIM_SUB: u64 = 29;
const PRIM_MUL: u64 = 31;
const PRIM_DIV: u64 = 33;
const PRIM_EQ: u64 = 35;
const PRIM_LT: u64 = 37;
const PRIM_GT: u64 = 39;
const PRIM_QUOTE: u64 = 41;
const PRIM_APP1: u64 = 43;

fn reverse_list(lst: Word) !Word {
    if (!is_linear(lst)) {
        var result: Word = NIL;
        var cur = lst;

        while (cur != NIL) {
            if (!is_list(cur)) {
                discard(result);
                return error.TypeMismatch;
            }
            const obj = get_ptr(cur);
            const elem_copy = try deep_copy(obj[0]);
            errdefer discard(elem_copy);
            const new_result = try make_list(elem_copy, result, true);
            errdefer {
                discard(new_result);
                discard(result);
            }
            result = new_result;
            cur = obj[1];
        }

        return result;
    }

    var result: Word = NIL;
    var cur = lst;

    while (cur != NIL) {
        if (!is_list(cur)) {
            discard(result);
            return error.TypeMismatch;
        }
        const obj = get_ptr(cur);
        const elem = obj[0];
        const next = obj[1];

        obj[0] = 0;
        obj[1] = 0;
        free_obj(obj);

        result = try make_list(elem, result, true);
        cur = next;
    }

    return result;
}

inline fn execute_primitive(code: u64) anyerror!void {
    switch (code) {
        PRIM_DUP => {
            if (sp <= btm) return error.StackUnderflow;
            const a = stack[sp - 1];
            if (is_linear(a)) {
                const copied = try deep_copy(a);
                errdefer discard(copied);
                try push(copied);
            } else {
                try push(a);
            }
        },
        PRIM_SWAP => {
            if (sp < btm + 2) return error.StackUnderflow;
            const tmp = stack[sp - 1];
            stack[sp - 1] = stack[sp - 2];
            stack[sp - 2] = tmp;
        },
        PRIM_POP => {
            if (sp <= btm) return error.StackUnderflow;
            sp -= 1;
            discard(stack[sp]);
        },
        PRIM_DIP => {
            if (sp < btm + 2) return error.StackUnderflow;
            const prog = pop();
            defer discard(prog);
            const x = pop();
            errdefer discard(x);

            try eval(prog);
            try push(x);
        },
        PRIM_ROT => {
            if (sp < btm + 3) return error.StackUnderflow;
            const a = stack[sp - 1];
            const b = stack[sp - 2];
            const c = stack[sp - 3];
            stack[sp - 3] = b;
            stack[sp - 2] = a;
            stack[sp - 1] = c;
        },
        PRIM_CONS => {
            if (sp < btm + 2) return error.StackUnderflow;
            const car_val = pop();
            errdefer discard(car_val);
            const cdr_val = pop();
            errdefer discard(cdr_val);
            const new_list = try make_list(car_val, cdr_val, true);
            push_unchecked(new_list);
        },
        PRIM_UNCONS => {
            if (sp <= btm) return error.StackUnderflow;
            const lst = pop();
            if (!is_list(lst)) {
                discard(lst);
                return error.TypeMismatch;
            }
            const obj = get_ptr(lst);
            const car_val = obj[0];
            errdefer discard(car_val);
            const cdr_val = obj[1];
            errdefer discard(cdr_val);

            if (is_linear(lst)) {
                obj[0] = 0;
                obj[1] = 0;
                free_obj(obj);
            }

            push_unchecked(cdr_val);
            try push(car_val);
        },
        PRIM_CAT => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = pop();
            errdefer discard(b);
            const a = pop();
            errdefer discard(a);

            if (!is_list(a) or !is_list(b)) {
                discard(a);
                discard(b);
                return error.TypeMismatch;
            }

            const a_rev = try reverse_list(a);
            errdefer discard(a_rev);

            var result = b;
            var cur = a_rev;

            while (cur != NIL) {
                if (!is_list(cur)) {
                    discard(cur);
                    discard(result);
                    return error.TypeMismatch;
                }
                const obj = get_ptr(cur);
                const elem = obj[0];
                const next = obj[1];

                obj[0] = 0;
                obj[1] = 0;
                free_obj(obj);

                result = try make_list(elem, result, true);
                cur = next;
            }

            push_unchecked(result);
        },
        PRIM_I => {
            if (sp <= btm) return error.StackUnderflow;
            const prog = pop();
            if (!is_list(prog)) {
                discard(prog);
                return error.TypeMismatch;
            }
            try eval(prog);
            discard(prog);
        },
        PRIM_FIRST => {
            if (sp <= btm) return error.StackUnderflow;
            const lst = pop();
            if (!is_list(lst)) {
                discard(lst);
                return error.TypeMismatch;
            }
            const obj = get_ptr(lst);
            const car_val = obj[0];

            if (is_linear(lst)) {
                obj[0] = 0;
                discard(obj[1]);
                obj[1] = 0;
                free_obj(obj);
            }

            push_unchecked(car_val);
        },
        PRIM_REST => {
            if (sp <= btm) return error.StackUnderflow;
            const lst = pop();
            if (!is_list(lst)) {
                discard(lst);
                return error.TypeMismatch;
            }
            const obj = get_ptr(lst);
            const cdr_val = obj[1];

            if (is_linear(lst)) {
                discard(obj[0]);
                obj[0] = 0;
                obj[1] = 0;
                free_obj(obj);
            }

            push_unchecked(cdr_val);
        },
        PRIM_IFTE => {
            if (sp < btm + 3) return error.StackUnderflow;
            const f_branch = pop();
            defer discard(f_branch);
            const t_branch = pop();
            defer discard(t_branch);
            const cond = pop();
            defer discard(cond);

            if (cond == 1) {
                try eval(f_branch);
            } else {
                try eval(t_branch);
            }
        },
        PRIM_BRANCH => {
            if (sp < btm + 3) return error.StackUnderflow;
            const f_branch = pop();
            errdefer discard(f_branch);
            const t_branch = pop();
            errdefer discard(t_branch);
            const cond = pop();
            defer discard(cond);

            if (!is_int(cond)) {
                discard(t_branch);
                discard(f_branch);
                return error.TypeMismatch;
            }

            if (get_int(cond) != 0) {
                discard(f_branch);
                push_unchecked(t_branch);
            } else {
                discard(t_branch);
                push_unchecked(f_branch);
            }
        },
        PRIM_ADD => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = a + b - 1;
        },
        PRIM_SUB => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = a - b + 1;
        },
        PRIM_MUL => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = ((a >> 1) * (b >> 1) << 1) | 1;
        },
        PRIM_DIV => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = ((a >> 1) / (b >> 1) << 1) | 1;
        },
        PRIM_EQ => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = if (a == b) T else F;
        },
        PRIM_LT => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = if (a < b) T else F;
        },
        PRIM_GT => {
            if (sp < btm + 2) return error.StackUnderflow;
            const b = stack[sp - 1];
            const a = stack[sp - 2];
            if ((a & b & 1) == 0) {
                discard(pop());
                discard(pop());
                return error.TypeMismatch;
            }
            sp -= 1;
            stack[sp - 1] = if (a > b) T else F;
        },
        PRIM_QUOTE => {
            if (sp <= btm) return error.StackUnderflow;
            const x = pop();
            errdefer discard(x);
            const quoted = try make_list(x, NIL, true);
            push_unchecked(quoted);
        },
        PRIM_APP1 => {
            if (sp < btm + 2) return error.StackUnderflow;
            const prog = pop();
            defer discard(prog);
            if (!is_list(prog)) {
                return error.TypeMismatch;
            }
            const btm_bk = btm;
            btm = sp - 1;
            defer btm = btm_bk;
            try eval(prog);
            if (btm != sp - 1) {
                return error.TypeMismatch;
            }
        },
        else => return error.UnknownPrimitive,
    }
}

fn init_symbols(a: Allocator) !void {
    sym_table = std.StringHashMap(Word).init(a);
    interned_strings = std.StringHashMap(Word).init(a);

    const syms = [_]struct { []const u8, u64 }{
        .{ "dup", PRIM_DUP },
        .{ "swap", PRIM_SWAP },
        .{ "pop", PRIM_POP },
        .{ "dip", PRIM_DIP },
        .{ "rot", PRIM_ROT },
        .{ "cons", PRIM_CONS },
        .{ "uncons", PRIM_UNCONS },
        .{ "cat", PRIM_CAT },
        .{ "i", PRIM_I },
        .{ "first", PRIM_FIRST },
        .{ "rest", PRIM_REST },
        .{ "ifte", PRIM_IFTE },
        .{ "branch", PRIM_BRANCH },
        .{ "+", PRIM_ADD },
        .{ "-", PRIM_SUB },
        .{ "*", PRIM_MUL },
        .{ "/", PRIM_DIV },
        .{ "=", PRIM_EQ },
        .{ "<", PRIM_LT },
        .{ ">", PRIM_GT },
        .{ "quote", PRIM_QUOTE },
        .{ "app1", PRIM_APP1 },
    };

    for (syms) |s| {
        const name = try intern_string(s[0]);
        const sym = try make_symbol(name, @as(Word, s[1]), false);
        try sym_table.put(s[0], sym);
    }
}

fn deinit_symbols() void {
    var it = interned_strings.iterator();
    while (it.next()) |entry| {
        gpa.free(entry.key_ptr.*);
    }
    interned_strings.deinit();
    sym_table.deinit();
}

// Evaluator
inline fn eval(prog: Word) anyerror!void {
    if (prog != NIL) {
        try eval_impl(prog);
    }
}

fn eval_impl(prog: Word) anyerror!void {
    @setEvalBranchQuota(10000);
    var cursor = prog;
    while (true) {
        if (sp >= STACK_MAX) return error.StackOverflow;
        inline for (0..9) |_| {
            const obj = get_ptr(cursor);
            const op = obj[0];
            cursor = obj[1];

            if (is_symbol(op)) {
                const op_val = symbol_op(op);
                if (is_int(op_val)) {
                    const code = @as(u64, op_val);
                    try execute_primitive(code);
                } else {
                    try eval(op_val);
                }
            } else if (is_linear(op)) {
                const copied = try deep_copy(op);
                try push(copied);
            } else {
                try push(op);
            }

            if (cursor == NIL) {
                if (sp >= STACK_MAX) return error.StackOverflow;
                return;
            }
        }
    }
}

inline fn eval_repl(op: Word) anyerror!void {
    if (sp >= STACK_MAX) return error.StackOverflow;
    switch (op & 0b1111) {
        TAG_SYM, TAG_SYM | FLAG_LINEAR => {
            const op_val = symbol_op(op);
            if (is_int(op_val)) {
                const code = @as(u64, op_val);
                try execute_primitive(code);
            } else {
                try eval(op_val);
            }
        },
        else => {
            try push(op);
        },
    }
}

// Parser
const Parser = struct {
    src: []const u8,
    pos: usize,
    make_persistent: bool,

    fn init(src: []const u8, make_persistent: bool) Parser {
        return Parser{ .src = src, .pos = 0, .make_persistent = make_persistent };
    }

    fn skip_whitespace(self: *Parser) void {
        while (self.pos < self.src.len and std.ascii.isWhitespace(self.src[self.pos])) {
            self.pos += 1;
        }
    }

    fn peek_token(self: *Parser) ?[]const u8 {
        const saved_pos = self.pos;
        defer self.pos = saved_pos;

        self.skip_whitespace();
        if (self.pos >= self.src.len) return null;

        const ch = self.src[self.pos];
        if (std.ascii.isWhitespace(ch) or ch == '[' or ch == ']') return null;

        const start = self.pos;
        while (self.pos < self.src.len and !std.ascii.isWhitespace(self.src[self.pos]) and
            self.src[self.pos] != '[' and self.src[self.pos] != ']')
        {
            self.pos += 1;
        }
        return self.src[start..self.pos];
    }

    fn parse_one(self: *Parser) anyerror!?Word {
        self.skip_whitespace();
        if (self.pos >= self.src.len) return null;

        if (self.peek_token()) |token| {
            if (std.mem.eql(u8, token, "DEFINE")) {
                try self.parse_define();
                return try self.parse_one();
            }
        }

        const ch = self.src[self.pos];

        if (ch == '[') {
            self.pos += 1;
            var items = std.ArrayList(Word).init(gpa);
            defer items.deinit();

            while (true) {
                self.skip_whitespace();
                if (self.pos >= self.src.len) return error.UnmatchedBracket;
                if (self.src[self.pos] == ']') {
                    self.pos += 1;
                    break;
                }
                const item = try self.parse_one() orelse return error.UnmatchedBracket;
                try items.append(item);
            }

            var result: Word = NIL;
            var i = items.items.len;
            while (i > 0) {
                i -= 1;
                result = try make_list(items.items[i], result, self.make_persistent == false);
            }
            return result;
        } else if (ch == ']') {
            return error.UnmatchedBracket;
        } else if (std.ascii.isDigit(ch) or (ch == '-' and self.pos + 1 < self.src.len and std.ascii.isDigit(self.src[self.pos + 1]))) {
            const start = self.pos;
            if (ch == '-') self.pos += 1;
            while (self.pos < self.src.len and std.ascii.isDigit(self.src[self.pos])) {
                self.pos += 1;
            }
            const num = try std.fmt.parseInt(i64, self.src[start..self.pos], 10);
            return make_int(num);
        } else {
            const start = self.pos;
            while (self.pos < self.src.len and !std.ascii.isWhitespace(self.src[self.pos]) and self.src[self.pos] != '[' and self.src[self.pos] != ']') {
                self.pos += 1;
            }
            const token = self.src[start..self.pos];
            const sym = lookup_symbol(token) orelse return error.UnknownSymbol;
            return sym;
        }
    }

    fn parse_define(self: *Parser) anyerror!void {
        self.skip_whitespace();
        while (self.pos < self.src.len and !std.ascii.isWhitespace(self.src[self.pos])) {
            self.pos += 1;
        }

        self.skip_whitespace();
        const name_start = self.pos;
        while (self.pos < self.src.len and !std.ascii.isWhitespace(self.src[self.pos]) and
            self.src[self.pos] != '[' and self.src[self.pos] != ']')
        {
            self.pos += 1;
        }
        if (name_start == self.pos) return error.InvalidDefine;
        const name = self.src[name_start..self.pos];

        self.skip_whitespace();
        if (self.pos + 1 >= self.src.len or
            self.src[self.pos] != '=' or self.src[self.pos + 1] != '=')
        {
            return error.InvalidDefine;
        }
        self.pos += 2;

        const name_word = try intern_string(name);
        const placeholder_sym = try make_symbol(name_word, NIL, false);
        try sym_table.put(name, placeholder_sym);

        var body = std.ArrayList(Word).init(gpa);
        defer body.deinit();

        var body_parser = Parser.init(self.src, true);
        body_parser.pos = self.pos;

        while (true) {
            body_parser.skip_whitespace();
            if (body_parser.pos >= body_parser.src.len) return error.InvalidDefine;

            if (body_parser.src[body_parser.pos] == '.') {
                body_parser.pos += 1;
                self.pos = body_parser.pos;
                break;
            }

            const item = try body_parser.parse_one() orelse return error.InvalidDefine;
            try body.append(item);
        }

        var result: Word = NIL;
        var i = body.items.len;
        while (i > 0) {
            i -= 1;
            result = try make_list(body.items[i], result, false);
        }

        const sym_obj = get_ptr(placeholder_sym);
        sym_obj[1] = result;
    }
};

fn print_word(w: Word) !void {
    if (is_int(w)) {
        try stdout.print("{d}", .{get_int(w)});
    } else if (is_list(w)) {
        try stdout.print("[", .{});
        var cur = w;
        var first = true;
        while (cur != NIL) {
            if (!is_list(cur)) break;
            if (!first) try stdout.print(" ", .{});
            first = false;
            try print_word(car(cur));
            cur = cdr(cur);
        }
        try stdout.print("]", .{});
    } else if (is_symbol(w)) {
        const obj = get_ptr(w);
        const name_word = obj[0];

        if (name_word != 0) {
            const name_ptr: [*:0]const u8 = @ptrFromInt(name_word);
            const name_slice = std.mem.span(name_ptr);
            try stdout.print("<sym: {s}>", .{name_slice});
        } else {
            try stdout.print("<sym>", .{});
        }
    } else {
        try stdout.print("[]", .{});
    }
}

fn print_stack() !void {
    try stdout.print("Stack: ", .{});
    for (stack[0..sp]) |w| {
        try print_word(w);
        try stdout.print(" ", .{});
    }
    try stdout.print("\n", .{});
}

pub fn main() !void {
    gpa = std.heap.page_allocator;
    try init_mem_pool();
    defer deinit_mem_pool();

    try init_symbols(gpa);
    defer deinit_symbols();

    defer cleanup_stack();

    const args = try std.process.argsAlloc(gpa);
    defer std.process.argsFree(gpa, args);

    var src: []const u8 = undefined;
    var should_free = false;
    var is_file = false;

    if (args.len > 1) {
        const file = try std.fs.cwd().openFile(args[1], .{});
        defer file.close();
        src = try file.readToEndAlloc(gpa, 1024 * 1024);
        should_free = true;
        is_file = true;
    } else {
        try stdout.print("joy-- > ", .{});
        src = try std.io.getStdIn().reader().readUntilDelimiterAlloc(gpa, '\n', 4096);
        should_free = true;
        is_file = false;
    }
    defer if (should_free) gpa.free(src);

    var parser = Parser.init(src, is_file);

    while (true) {
        const op = parser.parse_one() catch |err| {
            cleanup_stack();
            try stderr.print("Parse error: {}\n", .{err});
            return err;
        };

        if (op == null) break;

        eval_repl(op.?) catch |err| {
            cleanup_stack();
            try stderr.print("Evaluation error: {}\n", .{err});
            return err;
        };
    }

    try print_stack();
}
