module.exports = {
    branches: ['main'],  // 监控的分支，例如 "main"
    plugins: [
        '@semantic-release/commit-analyzer',       // 分析提交信息
        '@semantic-release/release-notes-generator', // 生成发布说明
        '@semantic-release/changelog',             // 更新 CHANGELOG.md 文件
        ['@semantic-release/npm', {
            npmPublish: false  // 如果要发布到 npm，请保持为 true
        }],
        ['@semantic-release/git', {
            assets: ['package.json', 'CHANGELOG.md'], // 提交更新的文件
            message: 'chore(release): ${nextRelease.version}' // 提交信息模板
        }],
        //'@semantic-release/github' // 在 GitHub 上发布 release（可选）
    ],
    preset: 'conventionalcommits',  // 使用 Conventional Commits 规范
    tagFormat: '${version}',        // 标签格式
};
